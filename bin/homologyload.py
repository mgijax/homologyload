#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create HomoloGene bcp files for the MRK_Cluster* tables
#
# Usage: homologyload.py
#
# Inputs:
#       1. Intermediate homologene file tab-delimited in following format:
#           1.
#           2.
#           3.
#           4.
#           5.
#           5.
#       2. Configuration - see homologyload.config
#         1.
#         2.
#         3.
#         4.
#         5.
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

TAB = '\t'
CRT = '\n'

# Lookup of members by hgID (cluster ID)
hgIdToMemberDict = {}

# Full paths to input and output files
inFile = os.environ['INPUT_FILE_LOAD']
clusterBCP = os.environ['CLUSTER_BCP']
memberBCP = os.environ['MEMBER_BCP']
accessionBCP = os.environ['ACCESSION_BCP']

fpInFile = ''
fpClusterBCP = ''
fpMemberBCP = ''
fpAccessionBCP = ''

nextClusterKey = ''
nextMemberKey = ''
nextAccessionKey = ''

clusterTypeKey = os.environ['CLUSTER_TYPE_KEY']
clusterSource = os.environ['CLUSTER_SRC_KEY']
clusterVersion =  os.environ['HOMOLOGY_VERSION']
clusterDate = os.environ['HOMOLOGY_DATE']

createdByKey = 1527
cdate = mgi_utils.date("%m/%d/%Y")

# for hgId to cluster Accession
ldbKey = os.environ['HOM_LDB_KEY']
mgiTypeKey =  os.environ['CLUSTER_MGITYPE_KEY']

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

    print 'Deleting Homology Cluster Members'
    db.sql('''select _Cluster_key
    into #todelete2
    from MRK_Cluster
    where _CreatedBy_key = %s''' % createdByKey, None)

    db.sql('create index idx1 on #todelete2(_Cluster_key)', None)

    db.sql('''delete MRK_ClusterMember
        from #todelete2 d, MRK_ClusterMember m
        where d._Cluster_key = m._Cluster_key''', None)

    db.sql('''delete MRK_Cluster
        from #todelete2 d, MRK_Cluster m
        where d._Cluster_key = m._Cluster_key''', None)

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

#############
# Main      #
#############
deleteHomologies()

print  'loading a dictionary mapping each hgID to its list of members'
for line in fpInFile.readlines():
    tokens = string.split(line[:-1], TAB)
    hgId = tokens[0]
    tokens = tokens[1:]
    if not hgIdToMemberDict.has_key(hgId):
	hgIdToMemberDict[hgId] = []   
    hgIdToMemberDict[hgId].append(tokens)

fpInFile.close()

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
    sequenceNum = 0
    for line in lineList:
	#print 'line: %s' % line
	sequenceNum += 1
	markerKey = line[2]
	#print 'markerKey: %s' % markerKey
	fpMemberBCP.write('%s%s%s%s%s%s%s%s' % (nextMemberKey, TAB, nextClusterKey, TAB, markerKey, TAB, sequenceNum, CRT))
	nextMemberKey += 1

    # now increment the cluster key
    nextClusterKey += 1

db.useOneConnection(0)
fpClusterBCP.close()
fpMemberBCP.close()
fpAccessionBCP.close()
