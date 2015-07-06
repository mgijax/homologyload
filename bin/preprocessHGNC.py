#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       From HGNC input file create load ready file
#
# Usage: preprocessHGNC.py
#
# Inputs:
#	1. HGNC file tab-delimited in following format:
#	    1. Human EntrezGene ID
#	    2. Mouse MGI IDs, comma-delimited
#           3. HGNC ID (sans the 'HGNC:')
#
#	2. Configuration - see hgncload.config
#
# Outputs:
#	 1. load ready file
#	 2. QC report file
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
import string
import mgi_utils
import clusterize
import db

db.setAutoTranslate()
db.setAutoTranslateBE()

###--- globals ---###

# constants
TAB= '\t'
CRT = '\n'

# EG ID/Human Marker associations from the database
# {egID:marker key, ...}
egToMarkerDict = {}

# MGI ID/Mouse Marker associations from the database
# {mgiID:marker key, ...}
mgiToMarkerDict = {}

#
# paths to input and output files
#

# input file from HGNC
inFilePath = os.environ['INPUT_FILE']

# path to the file to be clustered by the clusterizer
# this was for debugging this first clustered load. I chose to leave this code
# in for future debugging purposes
clustererFilePath = os.environ['INPUT_FILE_CLUSTERER']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

#
# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'
head = 'Line#%sHuman EG ID%sMouse MGI IDs%s HGNC ID%s' % (TAB, TAB, TAB, CRT)

rptOne = 'Lines where Human EG ID not in database%s%s%s%s' % (CRT, CRT, head, sep)

rptTwo = '%s%sLines where a Mouse MGI ID not in database%s%s%s%s' % (CRT, CRT,CRT, CRT, head, sep)

#
# file descriptors
#

fpInFile = ''

# idList replaces this file, but keeping for debugging purposes
fpClustererFile = ''
fpLoadFile = ''
fpQcRpt = ''

###--- functions ---###

def init():
    # Purpose: Initialization of  database connection and file descriptors,
    #       create database lookup dictionaries; create dictionary from
    #       input file
    # Returns: 1 if file descriptors cannot be initialized
    # Assumes: Nothing
    # Effects: opens a database connection
    # Throws: Nothing


    global egToMarkerDict, mgiToMarkerDict
    global fpInFile, fpClustererFile, fpLoadFile, fpQcRpt

    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    try:
	fpInFile = open(inFilePath, 'r')
    except:
	exit('Could not open file for reading %s\n' % inFilePath)
    # idList replaces this file, but keeping for time being
    try:
	fpClustererFile = open(clustererFilePath, 'w')
    except:
	exit('Could not open file for writing %s\n' % clustererFilePath)
    try:
        fpLoadFile = open(loadFilePath, 'w')
    except:
        exit('Could not open file for writing %s\n' % loadFilePath)

    try:
	fpQcRpt = open(qcRptPath, 'w')
    except:
	exit('Could not open file for writing %s\n' % qcRptPath)


    # get all human markers that are associated with egIds
    results = db.sql('''select distinct a.accid as egId, m._Marker_key
	from ACC_Accession a, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 55
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Status_key in (1,3)
	and m._Organism_key = 2''', 'auto')
    #
    # create Human egID to marker lookup from database
    #
    for r in results:
	egId = r['egId']
	markerKey = r['_Marker_key']
	egToMarkerDict[egId] = markerKey

    # get all mouse markers
    results = db.sql('''select distinct a.accID as mgiId, m._Marker_key
	    from ACC_Accession a, MRK_Marker m
	    where a._MGIType_key = 2
	    and a._LogicalDB_key = 1
	    and a.prefixPart = 'MGI:'
	    and a._Object_key = m._Marker_key
	    and m._Marker_Status_key in (1,3)
	    and m._Organism_key = 1''', 'auto')

    # removed per Richard
    # and a.preferred = 1
    for r in results:
	    mgiId = r['mgiId']
	    markerKey = r['_Marker_key']
	    mgiToMarkerDict[mgiId] = markerKey

    return

def process():
    # Purpose: Create load ready file from HGNC file and the database
    # Returns: 0
    # Assumes: All lookup structures have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    global rptOne, rptTwo, hgncIdOnlyCount, mgiIdOnlyCount
    
    # dictionary of id pairs to send to the clusterizer
    toClusterList = []
    lineCt = 0
    # These will be ignored by the load
    hgncIdOnlyCount = 0

    # These will be ignored by the load until we determine how to 
    # cluster them in the absence of mouse MGI IDs
    mgiIdOnlyCount = 0

    for line in fpInFile.readlines():
	lineCt += 1
	
	(egId, mgiIDs, junk1) = string.split(line[:-1], TAB)
	egId = string.strip(egId)

	# if both egId and mgiId columns are blank, skip and don't
	# report
	if egId == '' and mgiIDs == '':
	    hgncIdOnlyCount += 1
	    continue

	# if mgiId, and no egID skip and don't report
	elif egId == '':
	    mgiIdOnlyCount += 1
	    continue
	elif mgiIDs == '':
	    mgiIDs = 'None'

	# 1 means error on this line
	error = 0

	# idList replaces this file, but keeping for time being
	clusterFileLine = ''

	# current cluster - if there are no errors it will be added to
        # 'toClusterList'
        currentClusterList = []

	# report and skip lines with egId not in the database
	if egId and egId not in egToMarkerDict.keys():
	    rptOne = '%s%s%s%s' % (rptOne, lineCt, TAB, line)
	    # if egId not in database continue to next input line
	    continue
	mgiIDList = map(string.strip, string.split(mgiIDs, ','))

	for id in mgiIDList:
	    id = string.strip(id)
	    # report and skip lines with mgiId not in the database
	    if id != 'None' and id not in mgiToMarkerDict.keys(): 
		error = 1
		rptTwo = '%s%s%s%s' % (rptTwo, lineCt, TAB, line)
		# No need to check any more ids, get out of the loop
		break
	    else:
		currentClusterList.append([egId, id])
	    # idList replaces this file, but keeping for time being
	    clusterFileLine = ('%s%s%s%s%s' % \
		(clusterFileLine, egId, TAB, id, CRT))
	# if any mgi IDs not in database continue to next input line
	if error == 1:
		continue

	# no errors so append the next cluster
	toClusterList = toClusterList + currentClusterList

	# idList replaces this file, but keeping for time being
	fpClustererFile.write(clusterFileLine)
	# if we get here, we the egId is in the database and ALL the
	# mgiIds are in the database
		
    # idList replaces this file, but keeping for debugging purposes
    fpClustererFile.close()

    # clusterDict = clusterize.cluster(clustererFilePath, 'HGNC')
    clusterDict = clusterize.cluster(toClusterList, 'HGNC')

    # now resolve the ids to database keys; human and mouse gene keys
    for clusterId in clusterDict.keys():
	idTuple = clusterDict[clusterId]
	humanKeyList = []
	mouseKeyList = []
	for id in idTuple:
	    if id.startswith('MGI:'):
		mouseKeyList.append(str(mgiToMarkerDict[id]))
	    else:
		humanKeyList.append(str(egToMarkerDict[id]))
	# we want human before mouse for cluster member sequence numbering
	keyList = humanKeyList + mouseKeyList
	keyString = ', '.join(keyList)
	fpLoadFile.write('%s%s%s%s' % (clusterId, TAB, keyString, CRT))

    return

def writeReports():
    # Purpose: writes out all sections of the QC report
    # Returns: 0
    # Assumes: rptOneand rptTwo have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)

    return

def closeFiles():
    # Purpose: closes file descriptors and database connection
    # Returns: 0
    # Assumes: file descriptors have been initialized
    # Effects:  None
    # Throws: Nothing

    fpInFile.close()
    fpLoadFile.close()
    fpQcRpt.close()

    # close the database connection
    db.useOneConnection(0)
    
    return

###--- main program ---###

print '%s' % mgi_utils.date()

print 'initializing'
init()

print 'processing clusters'
process()

print 'writing reports'
writeReports()

print 'closing files'
closeFiles()

print '%s' % mgi_utils.date()
