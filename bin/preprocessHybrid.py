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
#	None
#
# Outputs:
#       List of connected components 
#
# Exit Codes:
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes:  None
#
#	Hybrid Cluster conflict values:
#	1. conflict
#	2. reject
#	3. none
#
#	Hybrid Cluster source values:
#	1. HNGC
#	2. HG
#	3. both
#	4. none
#
###########################################################################
import db
import string
import Set
import os

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
connCompList = []
connCompDict = {}
# List of the hybrid clusters
hybridClusterList = []

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
    global fpCC, fpH
    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)
    fpCC = open('connComp.rpt', 'w')
    fpH = open('hybridCluster.rpt', 'w')
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
	clusterDict[cKey].source = 'hgnc'
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
    print '# HGNC Clusters: %s' % len(hgncDict)
    print '# Markers found: %s' % len(keyToMarkerDict)


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
	clusterDict[cKey].source = 'homologene'
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
    print '# HG Clusters: %s' % len(homologeneDict)
    print '# Markers found: %s' % len(keyToMarkerDict)

    #
    # Create list of all clusters from both providers
    allClustersList = hgncDict.values() + homologeneDict.values()
    print 'Total HGNC/HG clusters found: %s' % len(allClustersList)

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
    #print 'closure source: %s' % source
    for m in c.markers:
	mKey = m.key
 	if source == 'homologene' and hgncDict.has_key(mKey):
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
	#print '\n' + c.toString()
	if c.done == 0:
	    connectedComp = []
            closure(c)
	    #connCompList.append(connectedComp) # replaced by dict below
	    hcList = decideWhatToDo(connectedComp)
	    hybridClusterList = hybridClusterList + hcList
	    # add to dict
            ccCount += 1
	    connCompDict[ccCount] = [connectedComp, hcList]

def writeComponents():
    global fpCC

    #ccCount = 0
    #for cc in connCompList:
	#ccCount += 1
    keyList = connCompDict.keys()
    keyList.sort()
    for ccCount in keyList:
	connCompList = connCompDict[ccCount][0]
	fpCC.write('comp%s:%s' % (ccCount, CRT))
	for c in connCompList:
	    fpCC.write('%s%s' % (c.toString(), CRT))
    fpCC.close()

def writeHybrid():
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
 
def decideWhatToDo(cc): # c is a connected componente; a list of Cluster objects

    #       Hybrid Cluster conflict values:
    #       1. conflict
    #       2. reject
    #       3. none
    #
    #       Hybrid Cluster source values:
    #       1. HNGC
    #       2. HG
    #       3. both
    #       4. none
    
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

    # Rule 3b only sgl species clusters, keep all HG, source=HG
    single = 1
    hgList = []

    for c in cc:
	if len(c.organisms) > 1:
	    single = 0
	    break
	if c.source == 'homologene':
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
	if c.source == 'homologene':
	    hgList.append(c)
	else:
	    hgncList.append(c)
	if len(c.organisms) == 2:
	    if c.source == 'homologene':
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
    #print 'bothHG: %s bothHGNC: %s' % (bothHG, bothHGNC)
    if bothHG == 1 and  bothHGNC == 1: 
	# Rule 1 - only one M/H cluster from each source; keep one source=both,
	# conflict=none
	#print 'len(hgList): %s len(hgncList): %s' % (len(hgList), len(hgncList))
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
	    #print CRT
	    #print 'hg: %s' % hgMarkerSet
	    #print 'hgnc: %s' % hgncMarkerSet
	    #print CRT
	    if not len(hgMarkerSet.difference(hgncMarkerSet)):
		#print 'rule 1 sets are equal'
		# both the same, arbitrarily pick hg

		# update the default hybridSource, conflict and rule values
		hgCluster.hybridSource = 'both'
		hgCluster.conflict = 'none'
		hgCluster.rule = '1'
		return [hgCluster]
	    else:
		# Rule 4 sources disagree keep hgnc cluster
		#print 'rule 4 one cluster each'
		hgncCluster.rule = '4'
		# hybridSource and conflict value remain the default
		return [hgncCluster]
	# Rule 4 both sources have disagreeing M/H cluster, keep hgnc clusters
	# source=hgnc, conflict= conflict
	else:
	    #print 'rule 4 multi clusters'
	    # update rule number from default to 4
	    for c in hgncList:
		c.rule = '4'
	    # hybridSource and conflict value remain the default
	    return hgncList
    #print "we shouldn't get here"
    return
#####################################
#
# Main
#
#####################################

init()
getClusters()
findComponents()
writeComponents()
writeHybrid()
