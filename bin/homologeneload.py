#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create HomoloGene (HG) bcp files for the MRK_Cluster* tables
#
# Usage: homologyload.py
#
# Inputs:
#       1. load-ready HG file tab-delimited in following format:
#           1. Cluster ID - identifies the HG cluster
#           2. Taxonomy ID - identifies the organism for this marker
#           3. EntrezGene ID - identifies the marker itself
#           4. Marker symbol
#           5. MGI ID - NCBI identifier for the RefSeq
#           6. Protein sequence ID - RefSeq
#       2. Configuration - see homologyload.config
#           1. INPUT_FILE_DEFAULT - Full path name HG file we copy to 
#		the load input directory
#	    2. INPUT_FILE_LOAD - Full path name of the load-ready file that is
#		created by this script
#	    3. CLUSTER_BCP - MRK_Cluster bcp file
#	    4. MEMBER_BCP - MRK_ClusterMember bcp file
#	    5. ACCESSION_BCP - ACC_Accession bcp file  ( HG ID assoc to cluster)
#	    6. CLUSTER_TYPE_KEY - MRK_Cluster._ClusterType_key
#	    7. CLUSTER_SRC_KEY - MRK_Cluster._ClusterSource_key
#	    8. RELEASE_NO_TEXT - Text to be prepended tO the HG release number
#	    9. HOMOLOGY_VERSION - HG release number stored with RELEASE_NO_TEXT
#		in MRK_Cluster.version
#	    10. HOM_LDB_KEY - HG logicalDB key for association HG ID with 
#		MRK_Cluster objects
#	    11. CLUSTER_MGITYPE_KEY - MGI Type key for MRK_Cluster
#	    12. MGD_DBUSER - database user loading this data
#	    13. MGD_DBPASSWORDFILE - pass work file containing pw for MGD_DBUSER
#
# Outputs:
#        1. MRK_Cluster.bcp
#        2. MRK_ClusterMember.bcp
#	 3. ACC_Accession.bcp
#
# Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Assumes:  Nothing
#
#  Notes:  None
#
###########################################################################
import sys
import os
import re
import db
import string
import mgi_utils
import loadlib
import symbolsort
import time

TAB = '\t'
CRT = '\n'

preferredOrganisms = [ 'human', 'mouse, laboratory', 'rat' ]

# Lookup of members by hgID (cluster ID)
hgIdToMemberDict = {}

# Full paths to input and output files
inFile = os.environ['INPUT_FILE_LOAD']
clusterBCP = os.environ['CLUSTER_BCP']
memberBCP = os.environ['MEMBER_BCP']
accessionBCP = os.environ['ACCESSION_BCP']

# file descriptors
fpInFile = ''
fpClusterBCP = ''
fpMemberBCP = ''
fpAccessionBCP = ''

# database primary keys, the next one available
nextClusterKey = ''	# MRK_Cluster
nextMemberKey = ''	# MRK_ClusterMember
nextAccessionKey = ''	# ACC_Accession

# get MRK_Cluster type and source keys from Configuration
clusterTypeKey = os.environ['CLUSTER_TYPE_KEY']
clusterSource = os.environ['CLUSTER_SRC_KEY']

# create MRK_Cluster.version for environment variabels
# RELEASE_NO_TEXT is set from config, HOMOLOGY_VERSION is set in calling script
clusterVersion = "%s %s" % (os.environ['RELEASE_NO_TEXT'], os.environ['HOMOLOGY_VERSION'])

# The timestamp on the HG input file is used for MRK_Cluster.cluster_date
clusterDate = time.strftime("%b %d, %Y",time.localtime(os.path.getmtime(os.environ['INPUT_FILE_DEFAULT'])))

# MGI_User key for this load
createdByKey = 1527
# today's date for record timestamp
cdate = mgi_utils.date("%m/%d/%Y")

# for HG ID to MRK_Cluster Accession
ldbKey = os.environ['HOM_LDB_KEY']
mgiTypeKey =  os.environ['CLUSTER_MGITYPE_KEY']

# create file descriptors for input/output files
try:
    fpInFile = open(inFile, 'r')
except:
    exit(1, 'Could not open file %s\n' % inFile)

try:
    fpClusterBCP = open(clusterBCP, 'w')
except:
    exit(1, 'Could not open file %s\n' % clusterBCP)

try:
    fpMemberBCP = open(memberBCP, 'w')
except:
    exit(1, 'Could not open file %s\n' % memberBCP)

try:
    fpAccessionBCP = open(accessionBCP, 'w')
except:
    exit(1, 'Could not open file %s\n' % accessionBCP)

# get next ACC_Accession, MRK_Cluster and MRK_ClusterMember key
user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
db.useOneConnection(1)
db.set_sqlUser(user)
db.set_sqlPasswordFromFile(passwordFileName)

results = db.sql('''select max(_Cluster_key) + 1 as nextKey
        from MRK_Cluster''', 'auto')
if results[0]['nextKey'] is None:
    nextClusterKey = 1000
else:
    nextClusterKey = results[0]['nextKey']

results = db.sql('''select max(_ClusterMember_key) + 1 as nextKey
        from MRK_ClusterMember''', 'auto')
if results[0]['nextKey'] is None:
    nextMemberKey = 1000
else:
    nextMemberKey = results[0]['nextKey']

results = db.sql('''select max(_Accession_key) + 1 as nextKey
        from ACC_Accession''', 'auto')
nextAccessionKey = results[0]['nextKey']

def deleteHomologies():
    # Purpose: delete accession, cluster and member records
    # Returns: nothing
    # Assumes: nothing
    # Effects: queries a database, deletes records from a database
    # Throws:  db.error, db.connection_exc

    print 'Deleting Homology Cluster Accessions'
    db.sql('''select _Accession_key
    into #todelete1
    from ACC_Accession
    where _CreatedBy_key = %s''' % createdByKey, None)

    db.sql('create index idx1 on #todelete1(_Accession_key)', None)

    db.sql('''delete ACC_Accession
        from #todelete1 d, ACC_Accession a
        where d._Accession_key = a._Accession_key''', None)

    db.sql('''select _Cluster_key
    into #todelete2
    from MRK_Cluster
    where _CreatedBy_key = %s''' % createdByKey, None)

    db.sql('create index idx1 on #todelete2(_Cluster_key)', None)
   
    print 'Deleting Homology Cluster Members' 
    db.sql('''delete MRK_ClusterMember
        from #todelete2 d, MRK_ClusterMember m
        where d._Cluster_key = m._Cluster_key''', None)

    print 'Deleting Homology Clusters'
    db.sql('''delete MRK_Cluster
        from #todelete2 d, MRK_Cluster m
        where d._Cluster_key = m._Cluster_key''', None)

def byOrganismAndSymbol (a, b):
     global x

     [aOrg, aSym] = a[:2]
     [bOrg, bSym] = b[:2]

     # easy cases first...

     # matching organisms, so purely a symbol sort
     if (aOrg == bOrg):
           return symbolsort.nomenCompare (aSym, bSym)

     # 'a' is from a preferred organism and 'b' is not
     if (aOrg in preferredOrganisms) and (bOrg not in preferredOrganisms):
           return -1

     # 'b' is from a preferred organism and 'a' is not
     if (aOrg not in preferredOrganisms) and (bOrg in preferredOrganisms):
           # 'b' comes first
           return 1

     # somewhat harder cases next...

     # both organisms are preferred (and they do not match), so need to
     # sort by position in the preferredOrganisms list

     if (aOrg in preferredOrganisms) and (bOrg in preferredOrganisms):
           aOrgPosition = preferredOrganisms.index(aOrg)
           bOrgPosition = preferredOrganisms.index(bOrg)
           return cmp(aOrgPosition, bOrgPosition)

     # neither organism is preferred (and they do not match), so need to
     # sort alphabetically by organism

     return cmp(aOrg, bOrg)

#############
# Main      #
#############
deleteHomologies()

# load input file into memory organizing by HG ID
for line in fpInFile.readlines():
    tokens = string.split(line[:-1], TAB)
    hgId = tokens[2]
    if not hgIdToMemberDict.has_key(hgId):
	hgIdToMemberDict[hgId] = []   
    hgIdToMemberDict[hgId].append(tokens)

fpInFile.close()

print "Creating BCP files"
for hgId in hgIdToMemberDict.keys():
    # lineList is list of lists [ [line1], [line2], ... ]
    lineList = hgIdToMemberDict[hgId]
    #    
    # create MRK_Cluster
    #
    fpClusterBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextClusterKey, TAB, clusterTypeKey, TAB, clusterSource, TAB, hgId, TAB, clusterVersion, TAB, clusterDate, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT))

    #
    # create ACC_Accession
    #
    fpAccessionBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextAccessionKey, TAB, hgId, TAB, TAB, hgId, TAB, ldbKey, TAB, nextClusterKey, TAB, mgiTypeKey, TAB, 0, TAB, 1, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT) )
    nextAccessionKey += 1

    #
    # create MRK_ClusterMember
    #

    # sort the list by organism and symbol
    lineList.sort(byOrganismAndSymbol)
    sequenceNum = 0
    for line in lineList:
	sequenceNum += 1
	markerKey = line[4]
	fpMemberBCP.write('%s%s%s%s%s%s%s%s' % (nextMemberKey, TAB, nextClusterKey, TAB, markerKey, TAB, sequenceNum, CRT))
	nextMemberKey += 1

    # now increment the cluster key
    nextClusterKey += 1

db.useOneConnection(0)
fpClusterBCP.close()
fpMemberBCP.close()
fpAccessionBCP.close()
