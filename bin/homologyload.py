
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
#       2. Configuration - see homologyload.config and individual load configs
#
# Outputs:
#        1. MRK_Cluster.bcp
#        2. MRK_ClusterMember.bcp
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
# sc   01/14/2015
#       - initial implementation
###########################################################################
import os
import mgi_utils
import time
import db

###--- globals ---###

# constants
TAB = '\t'
CRT = '\n'

# Full paths to input and output files
inFile = os.environ['INPUT_FILE_LOAD']
clusterBCP = os.environ['CLUSTER_BCP']
memberBCP = os.environ['MEMBER_BCP']

# file descriptors
fpInFile = ''
fpClusterBCP = ''
fpMemberBCP = ''

# database primary keys, the next one available
nextClusterKey = ''	# MRK_Cluster
nextMemberKey = ''	# MRK_ClusterMember

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

###--- functions ---###

def init():
    # Purpose: Initialization of  database connection and file descriptors,
    #       create database lookup dictionaries; create dictionary from
    #       input file
    # Returns: 1 if file descriptors cannot be initialized
    # Assumes: Nothing
    # Effects: opens a database connection
    # Throws: Nothing

    global fpInFile, fpClusterBCP, fpMemberBCP
    global createdByKey, nextClusterKey, nextMemberKey

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

    # get next MRK_Cluster and MRK_ClusterMember key
    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    results = db.sql('''select _User_key
            from MGI_User
            where login = '%s' ''' % createdBy, 'auto')

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

    return

def deleteHomologies():
    # Purpose: delete cluster and member records
    # Returns: nothing
    # Assumes: nothing
    # Effects: queries a database, deletes records from a database
    # Throws:  db.error, db.connection_exc

    db.sql('''select _Cluster_key
    into temp todelete2
    from MRK_Cluster
    where _CreatedBy_key = %s''' % createdByKey, None)

    db.sql('create index todelete2_idx2 on todelete2(_Cluster_key)', None)
   
    print('Deleting Homology Clusters, Members, Properties and Accessions') 
    db.sql('''delete from MRK_Cluster m
        using todelete2 d
        where d._Cluster_key = m._Cluster_key''', None)

    db.commit()

    return

def createBCPFiles():
    # Purpose: Create bcp files from the load ready file 
    # Returns: 0
    # Assumes: Nothing
    # Effects: Writes to the file system
    # Throws: Nothing

    global nextClusterKey,  nextMemberKey
    for line in fpInFile.readlines():
        tokens = str.split(line[:-1], TAB)
        id = tokens[0]
        members = tokens[1]
        memberList = list(map(str.strip, str.split(members, ',')))
        properties = ''
        if len(tokens) == 3:
            properties = tokens[2]

        #
        # create MRK_Cluster
        #
        fpClusterBCP.write('%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s' % (nextClusterKey, TAB, clusterTypeKey, TAB, clusterSource, TAB, TAB, TAB, clusterDate, TAB, createdByKey, TAB, createdByKey, TAB, cdate, TAB, cdate, CRT))

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

    return

def closeFiles():
    # Purpose: closes file descriptors and database connection
    # Returns: 0
    # Assumes: file descriptors have been initialized
    # Effects:  None
    # Throws: Nothing

    db.useOneConnection(0)
    fpClusterBCP.close()
    fpMemberBCP.close()
    return

###--- main program ---###

print('%s' % mgi_utils.date())

init()
deleteHomologies()
createBCPFiles()
closeFiles()

print('%s' % mgi_utils.date())
