#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create HomoloGene (HG) bcp files for the MRK_Cluster* tables
#
# Usage: homologeneload.py
#
# Inputs:
#       1. load-ready HG file tab-delimited in following format:
#           1. Cluster ID - identifies the HG cluster
#           2. Taxonomy ID - identifies the organism for this marker
#           3. EntrezGene ID - identifies the marker itself
#           4. Marker symbol
#           5. MGI ID - NCBI identifier for the RefSeq
#           6. Protein sequence ID - RefSeq
#       2. Configuration - see homologeneload.config
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
#  Notes: With the scrum-bob Misc project TR11886, we added 5 new homology 
#    loads to this product. This is the original script 
#    (formerly called homologyload.py) that ran the homologeneload. We 
#    renamed it and continue to use it to run the homologene load (see LOADER 
#    in the homologeneload.config) because it has many HomoloGene specific 
#    attributes.
#
#  We then wrote a generic	homologyload.py which all other homology 
#    loads use (see LOADER in individual load configs)
#
# History:
#
# sc   01/14/2015
#       - initial implementation
###########################################################################
import os
import string
import mgi_utils
import symbolsort
import time
import db

#db.setAutoTranslate(False)
#db.setAutoTranslateBE(False)

###--- globals ---###

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

def init():
    # Purpose: Initialization of  database connection and file descriptors,
    #       and next available database keys
    # Returns: 1 if file descriptors cannot be initialized
    # Assumes: Nothing
    # Effects: opens a database connection
    # Throws: Nothing

    global fpInFile, fpClusterBCP, fpMemberBCP, fpAccessionBCP 
    global nextClusterKey, nextMemberKey, nextAccessionKey

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

    return

def deleteHomologies():
    # Purpose: delete accession, cluster and member records
    # Returns: nothing
    # Assumes: nothing
    # Effects: queries a database, deletes records from a database
    # Throws:  db.error, db.connection_exc

    db.sql('''select _Cluster_key
    into #todelete2
    from MRK_Cluster
    where _CreatedBy_key = %s''' % createdByKey, None)

    db.sql('create index todelete2_idx1 on #todelete2(_Cluster_key)', None)
   
    print 'Deleting Homology Clusters and Members; any Properties and Accessions'
    db.sql('''delete from MRK_Cluster m
        using #todelete2 d
        where d._Cluster_key = m._Cluster_key''', None)

    db.commit()

    return

def byOrganismAndSymbol (a, b):
    # Purpose: function to pass to python list sort to sort 
    # cluster members so they get the proper sequence number
    # Returns: 1 if file descriptors cannot be initialized
    # Assumes: a and b have at least two elements
    # Effects: Nothing
    # Throws: Nothing

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

def process():
    # Purpose: Create bcp files from the load ready HomoloGene file  created
    #   by the preprocessor
    # Returns: 0
    # Assumes: All lookup structures have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    global hgIdToMemberDict, nextClusterKey, nextAccessionKey, nextMemberKey

    # load input file into memory organizing by HG ID
    for line in fpInFile.readlines():
	tokens = string.split(line[:-1], TAB)
	hgId = tokens[2]
	if not hgIdToMemberDict.has_key(hgId):
	    hgIdToMemberDict[hgId] = []
	hgIdToMemberDict[hgId].append(tokens)

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
	    fpMemberBCP.write('%s%s%s%s%s%s%s%s' % \
		(nextMemberKey, TAB, nextClusterKey, TAB, markerKey, TAB, sequenceNum, CRT))
	    nextMemberKey += 1

	# now increment the cluster key
	nextClusterKey += 1

    return

def closeFiles():

    fpInFile.close()
    fpClusterBCP.close()
    fpMemberBCP.close()
    fpAccessionBCP.close()
    db.useOneConnection(0)

    return

###--- main program ---###

init()
deleteHomologies()
process()

