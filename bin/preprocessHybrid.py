#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Find HGNC/HomoloGene connected components from MRK_Cluster tables
#	Apply rules to pick representative cluster
#	Create homology cluster input file	
#
# Usage: preprocessHybrid.py
#
# Inputs:
#	The MGD database
#
# Outputs:
#      1. load ready file
#      2. QC report file
#
# Exit Codes:
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes: 
#
#	Hybrid Cluster source property values:
#	1. HGNC
#	2. HomoloGene
#	3. HGNC and HomoloGene
#
#	Hybrid Cluster conflict values: 
#	  This script reflects an outdated conflict algorithm as it 
#	  was implemented, then decided we wouldn't load the conflict property
#
###########################################################################
import string
import Set
import os
import mgi_utils

try:
    if os.environ['DB_TYPE'] == 'postgres':
        import pg_db
        db = pg_db
        db.setTrace()
        db.setAutoTranslateBE()
    else:
        import db

except:
    import db

CRT = '\n'
TAB = '\t'

# HGNC clusters from the database {markerKey:Cluster, ...}
hgncDict = {}

# HomoloGene clusters from the database
homologeneDict = {}

# all cluster marker symbols by marker key
keyToMarkerDict = {}

# list of all Cluster objects from both providers
#
allClustersList = []

# list of all connected components
# {ccID: [ [the connected component], [the hybrid clusters]], ...}
connCompDict = {}

# List of the hybrid clusters
hybridClusterList = []

# load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

class Marker:
    # Is: data object  representing a Marker
    # Is: data object  representing a marker attributes
    # Does: provides direct access to its attributes
    #
    def __init__ (self):
        # Purpose: constructor
        # Returns: nothing
        # Assumes: nothing
        # Effects: nothing
        # Throws: nothing
	self.key = None
	self.symbol = None
	self.organism = None
    def toString(self):
	return '''%s|%s|%s''' % (self.key, self.symbol, self.organism)


class Cluster:
    # Is: data object  representing a homology cluster
    # Has: a set of cluster attributes
    # Does: provides direct access to its attributes
    #
    def __init__ (self):
        # Purpose: constructor
        # Returns: nothing
        # Assumes: nothing
        # Effects: nothing
        # Throws: nothing
        self.clusterKey = None

  	self.source = None

	self.done = 0
	self.markers = [] # list of marker objects
	self.organisms = set([]) # 'm' or 'h' 
	# for hybrid clusters
	self.conflict = None
	self.hybridSource = None
 	self.rule = None

    def toLoadFormat(self):
	humanKeyList = []
	mouseKeyList = []
	clusterList = []
	for m in self.markers:
	    if m.organism == 'h':
		humanKeyList.append('%s' % m.key)
	    else:
		mouseKeyList.append('%s' % m.key)
	# we want human before mouse for cluster member sequence numbering
	keyList = humanKeyList + mouseKeyList
	keyString = ', '.join(keyList)
	propertyString = 'secondary source:%s' % self.hybridSource
	clusterList.append('%s%s%s%s%s%s' % (self.clusterKey, TAB, keyString, TAB, propertyString, CRT))
	return ''.join(clusterList)

    def toString(self):
	mList = []
	for m in self.markers:
	    mList.append(m.toString())
	return 'cKey: %s%ssource: %s%smembers:%s%s%s' % (self.clusterKey, CRT, self.source, CRT, CRT, string.join(mList, CRT ), CRT)

    def toStringHybrid(self):
	mList = []
        for m in self.markers:
            mList.append(m.toString())
	return 'originalClustkey: %s%srule:%s%ssource: %s%shybrid source: %s%sconflict: %s%smembers:%s%s%s' % (self.clusterKey, CRT, self.rule, CRT, self.source, CRT, self.hybridSource, CRT, self.conflict, CRT, CRT, string.join(mList, CRT ), CRT)

def init():
    global fpCC, fpH, fpRptFile, fpLoadFile
    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
    # temp report
    fpCC = open('connComp.rpt', 'w')
    # temp report
    fpH = open('hybridCluster.rpt', 'w')
    # Sue's temp report 
    fpRptFile = open(os.environ['HYBRID_RPT'], 'w')
    # Load file in homologyload format
    try:
        fpLoadFile = open(loadFilePath, 'w')
    except:
        exit('Could not open file for writing %s\n' % loadFilePath)

    return 0

def getClusters():
    global hgncDict, homologeneDict, keyToMarkerDict, allClustersList

    #
    # Create HGNC Clusters
    #
   
    # HGNC is only human and mouse 
    results = db.sql('''select mcm._Cluster_key, mcm._Marker_key, m.symbol,
	    o.commonName
	from MRK_Cluster mc, MRK_ClusterMember mcm, MRK_Marker m, MGI_Organism o
	where mc._ClusterType_key = 9272150
	and mc._ClusterSource_key = 13437099
	and mc._Cluster_key = mcm._Cluster_key
	and mcm._Marker_key = m._Marker_key
	and m._Organism_key = o._Organism_key''', 'auto')

    # create HGNC cluster objects
    clusterDict = {}
    for r in results:
	cKey = r['_Cluster_key']
	mKey = r['_Marker_key']
	symbol = r['symbol']
	keyToMarkerDict[mKey] = symbol

	if not clusterDict.has_key(cKey):
	   clusterDict[cKey] = Cluster()
        clusterDict[cKey].clusterKey = cKey
	clusterDict[cKey].source = 'HGNC'
	m = Marker()
	m.key = mKey
	m.symbol = symbol
	m.organism = r['commonName']
	clusterDict[cKey].markers.append(m)
	o = 'm'
	if m.organism == 'human':
	    o = 'h'
        clusterDict[cKey].organisms.add(o)
	     

    # now map the markers in each cluster to its cluster
    for cKey in clusterDict.keys():
	currentCluster = clusterDict[cKey]
	for m in currentCluster.markers:
	    hgncDict[m.key] = currentCluster


    #
    # Create HomoloGene Clusters
    #

    # select only human and mouse
    results = db.sql('''select mcm._Cluster_key, mcm._Marker_key, m.symbol,
	    o.commonName
        from MRK_Cluster mc, MRK_ClusterMember mcm, MRK_Marker m, MGI_Organism o
        where mc._ClusterType_key = 9272150
        and mc._ClusterSource_key = 9272151
        and mc._Cluster_key = mcm._Cluster_key
        and mcm._Marker_key = m._Marker_key
	and m._Organism_key in (1,2) 
	and m._Organism_key = o._Organism_key''', 'auto')

   # create HomoloGene cluster objects
    clusterDict = {}
    for r in results:
        cKey = r['_Cluster_key']
        mKey = r['_Marker_key']
	symbol = r['symbol']
        keyToMarkerDict[mKey] = symbol

        if not clusterDict.has_key(cKey):
           clusterDict[cKey] = Cluster()
        clusterDict[cKey].clusterKey = cKey
	clusterDict[cKey].source = 'HomoloGene'
	m = Marker()
        m.key = mKey
        m.symbol = symbol
        m.organism = r['commonName']
        clusterDict[cKey].markers.append(m)
	o = 'm'
        if m.organism == 'human':
            o = 'h'
        clusterDict[cKey].organisms.add(o)

    # now map the markers in each cluster to its cluster
    for cKey in clusterDict.keys():
        currentCluster = clusterDict[cKey]
        for m in currentCluster.markers:
            homologeneDict[m.key] = currentCluster

    # Create list of all clusters from both providers
    allClustersList = hgncDict.values() + homologeneDict.values()

def closure(c): # c is Cluster object
    global connectedComp
    # do depth first search starting at cluster c
    # and set connectedComponent to the set of clusters in the connected 
    # component containing c
    if c.done == 1:
	return 0
    c.done = 1
    connectedComp.append(c)
    source = c.source
    for m in c.markers:
	mKey = m.key
 	if source == 'HomoloGene' and hgncDict.has_key(mKey):
	    c2 = hgncDict[mKey]
	    closure(c2)
	else:
	   if homologeneDict.has_key(mKey):
		c2 = homologeneDict[mKey]
		closure(c2)

def findComponents():
    global connectedComp, hybridClusterList, connCompDict
    ccCount = 0
    for c in allClustersList:
	if c.done == 0:
	    connectedComp = []
	    ccCount += 1
            closure(c)
	    #connCompList.append(connectedComp) # replaced by dict below
	    hcList = decideWhatToDo(connectedComp)
	    hybridClusterList = hybridClusterList + hcList

	    # add to dict
	    connCompDict[ccCount] = [connectedComp, hcList]
	    writeLoadFile(hcList)

def writeLoadFile(hcList): # h is list of Cluster objects from a given cc
    for c in hcList:
	fpLoadFile.write(c.toLoadFormat())

# write my temp cc report
def writeMyComponents():
    global fpCC

    keyList = connCompDict.keys()
    keyList.sort()
    for ccCount in keyList:
	connCompList = connCompDict[ccCount][0]
	fpCC.write('comp%s:%s' % (ccCount, CRT))
	for c in connCompList:
	    fpCC.write('%s%s' % (c.toString(), CRT))
    fpCC.close()

# write my temp hybrid report
def writeMyHybrid():
    global fpH

    #for c in hybridClusterList:
    keyList = connCompDict.keys()
    keyList.sort()
    for ccCount in keyList:
	hybridList = connCompDict[ccCount][1]
	fpH.write('comp%s:%s' % (ccCount, CRT))
	for c in hybridList:
	    fpH.write('%s%s' % (c.toStringHybrid(), CRT))
    fpH.close()


def writeHybrid():
    global fpRptFile

    keyList = connCompDict.keys()
    keyList.sort()
    for ccCount in keyList:
	ccID = 'comp%s' % ccCount
	ccClusterList = connCompDict[ccCount][0]
	hgList = []
	hgncList = []
	hybridSource = ''
	conflict = ''
	rule = ''
	compClustMarkerList = []
	hybridClustMarkerList = []
	for c in ccClusterList:
	    markerList = []
	    for m in c.markers:
		markerList.append(m.symbol)
	    markerString = '(%s)' % ', '.join(markerList)
	    compClustMarkerList.append(markerString)
	    if c.source == 'HomoloGene':
		hgList.append(markerString)
	    else:
		hgncList.append(markerString)
	hcList = connCompDict[ccCount][1]
	for hc in hcList:
	    hybridSource = hc.hybridSource
	    conflict = hc.conflict
	    rule = hc.rule
	    hybridList = []
	    for m in hc.markers:
		hybridList.append(m.symbol)
	    hString = '(%s)' % ', '.join(hybridList)
	    hybridClustMarkerList.append(hString)
	clusterString =  ', '.join(compClustMarkerList)
	hgncString = ', '.join(hgncList)
	hgString = ', '.join(hgList)
	hybridString = ', '.join(hybridClustMarkerList)
	fpRptFile.write('%s##%s##%s##%s##%s##%s##%s##%s%s' % (ccID, clusterString, hgncString, hgString, hybridString, hybridSource, conflict, rule, CRT))

def decideWhatToDo(cc): # c is a connected componente; a list of Cluster objects

    # Rule 2 - one cluster (so one source) --> keep it; conflict='none'
    if len(cc) == 1:
	c = cc[0]
	c.hybridSource = c.source
	c.conflict = 'none'
	c.rule = '2'
        return [c]
	

    #
    # multi sources
    #

    # Rule 3b only sgl species clusters, keep all HG, source=HomoloGene
    single = 1
    hgList = []

    for c in cc:
	if len(c.organisms) > 1:
	    single = 0
	    break
	if c.source == 'HomoloGene':
	    c.hybridSource = c.source
	    c.conflict = 'none'
	    c.rule = '3b'
	    hgList.append(c)
    if single == 1:	
	return hgList
    #
    # multi sources at least one multi species cluster
    #
    
    bothHGNC = 0
    bothHG = 0
    hgncList = []
    hgList = []
    for c in cc:
	# default to source=hybridSource, conflict= 'conflict'
	c.hybridSource = c.source
	c.conflict = 'conflict'
	c.rule = '3'
	if c.source == 'HomoloGene':
	    hgList.append(c)
	else:
	    hgncList.append(c)
	if len(c.organisms) == 2:
	    if c.source == 'HomoloGene':
		bothHG = 1
	    else:
		bothHGNC = 1
	
    # Rule 3 - Only one source has M/H clusters, keep clusters from that source
    # conflict='conflict
    if not (bothHG == bothHGNC): # if only one source has M/H
	# hybridSource and conflict values remain the default
	if bothHG:
	    return hgList
	else:
	    return hgncList

    # both sources have M/H clusters
    # don't need this test - it is implied, test w/o later
    if bothHG == 1 and  bothHGNC == 1: 
	# Rule 1 - only one M/H cluster from each source; keep one source=both,
	# conflict=none
	if len(hgList) == 1 and len(hgncList) == 1:
	    hgMarkerList = []
	    hgncMarkerList = []
	    hgCluster = hgList[0]
	    hgncCluster = hgncList[0]
	    for m in hgCluster.markers:
		hgMarkerList.append(m.key)
	    for m in hgncCluster.markers:
		hgncMarkerList.append(m.key)
	    hgMarkerSet = set(hgMarkerList)
	    hgncMarkerSet = set(hgncMarkerList)
	    if not len(hgMarkerSet.difference(hgncMarkerSet)):
		# both the same, arbitrarily pick hg
		# update the default hybridSource, conflict and rule values
		hgCluster.hybridSource = 'HomoloGene and HGNC'
		hgCluster.conflict = 'none'
		hgCluster.rule = '1'
		return [hgCluster]
	    else:
		# Rule 4 sources disagree keep hgnc cluster
		hgncCluster.rule = '4'
		# hybridSource and conflict value remain the default
		return [hgncCluster]
	# Rule 4 both sources have disagreeing M/H cluster, keep hgnc clusters
	# source=hgnc, conflict= conflict
	else:
	    # update rule number from default to 4
	    for c in hgncList:
		c.rule = '4'
	    # hybridSource and conflict value remain the default
	    return hgncList
    return

def closeFiles():
    fpLoadFile.close()
    fpRptFile.close()

    # close the database connection
    db.useOneConnection(0)

    return

#####################################
#
# Main
#
#####################################
print '%s' % mgi_utils.date()
print 'initializing'
init()

print 'getting clusters from the database'
getClusters()

print 'processing clusters'
findComponents()

print 'writing reports'
writeMyComponents()
writeMyHybrid()
writeHybrid()

print 'closing files'
closeFiles()

print '%s' % mgi_utils.date()

