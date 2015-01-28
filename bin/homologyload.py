#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create bcp files for the MRK_Cluster* tables
#
# Usage: homologyload.py
#
# Inputs:
#       1. load-ready file tab-delimited in following format:
#           1. Cluster ID - identifies the HG cluster
#           2. Comma separated list of marker keys in the cluster
#		ordered as you would like them sequenced in MGI_ClusterMember
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
#	    8. CLUSTER_MGITYPE_KEY - MGI Type key for MRK_Cluster
#	    9. MGD_DBUSER - database user loading this data
#	    10. MGD_DBPASSWORDFILE - pass work file containing pw for MGD_DBUSER
#	    11. LOAD_ACCESSION - true if loading cluster accIDs
# Outputs:
#        1. MRK_Cluster.bcp
#        2. MRK_ClusterMember.bcp
#	 3. ACC_Accession.bcp (optionally)
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
import time

TAB = '\t'
CRT = '\n'

# if true we are loading cluster accession IDs
loadAccessions = os.environ['LOAD_ACCESSION']

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

# The timestamp on the HG input file is used for MRK_Cluster.cluster_date
clusterDate = time.strftime("%b %d, %Y",time.localtime(os.path.getmtime(os.environ['INPUT_FILE_DEFAULT'])))

# MGI_User key for this load
createdBy = os.environ['JOBSTREAM']
createdByKey = ''

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

if loadAccessions == 'true':
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

results = db.sql('''select _User_key
	from MGI_User
	where login = "%s"''' % createdBy, 'auto')

createdByKey = results[0]['_User_key']
#print 'createdByKey: %s' % createdByKey
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

#############
# Main      #
#############
deleteHomologies()

print "Creating BCP files"
for line in fpInFile.readlines():
    id, members = string.split(line[:-1], TAB)
    memberList = string.split(members, ',')

    #
    # create MRK_Cluster
    #
    if loadAccessions == 'true':
	fpClusterBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextClusterKey, TAB, clusterTypeKey, TAB, clusterSource, TAB, id, TAB, TAB, clusterDate, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT))
    else:
	fpClusterBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextClusterKey, TAB, clusterTypeKey, TAB, clusterSource, TAB, TAB, TAB, clusterDate, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT))

    if loadAccessions == 'true':
	#
	# create ACC_Accession
	#
	fpAccessionBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextAccessionKey, TAB, id, TAB, TAB, id, TAB, ldbKey, TAB, nextClusterKey, TAB, mgiTypeKey, TAB, 0, TAB, 1, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT) )
	nextAccessionKey += 1

    #
    # create MRK_ClusterMember
    #

    # sort the list by organism and symbol
    sequenceNum = 0
    for markerKey in memberList:
        sequenceNum += 1
        fpMemberBCP.write('%s%s%s%s%s%s%s%s' % (nextMemberKey, TAB, nextClusterKey, TAB, markerKey, TAB, sequenceNum, CRT))
        nextMemberKey += 1

    # now increment the cluster key
    nextClusterKey += 1

db.useOneConnection(0)
fpClusterBCP.close()
fpMemberBCP.close()
if loadAccessions == 'true':
    fpAccessionBCP.close()
