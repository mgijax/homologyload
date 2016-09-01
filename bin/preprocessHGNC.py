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
db.setAutoTranslate(False)
db.setAutoTranslateBE(False)

###--- globals ---###

# constants
TAB= '\t'
CRT = '\n'

# Human EG ID/Mouse MGI ID associations from the file
# {egID:[list of mouse MGI IDs], ...}
humanEgToMouseMgiDict = {}

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
head = 'Line#%sHuman EG ID%sMouse MGI IDs%s' % (TAB, TAB, CRT)

rptOne = 'Lines where Human EG ID not in database (excludes skipped lines)%s%s%s%s' % (CRT, CRT, head, sep)

rptTwo = '%s%sLines where a Mouse MGI ID not in database (excludes skipped lines)%s%s%s%s' % (CRT, CRT,CRT, CRT, head, sep)
rptDebug = '%s%sInput resolved to keys%s%s' % (CRT, CRT, CRT, sep)
#
# file descriptors
#

fpInFile = ''

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


    # get all human markers that are associated with egIDs
    results = db.sql('''select distinct a.accid as egID, m._Marker_key
	from ACC_Accession a, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 55
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Status_key = 1
	and m._Organism_key = 2''', 'auto')
    #
    # create Human egID to marker lookup from database
    #
    for r in results:
	egID = r['egID']
	markerKey = r['_Marker_key']
	egToMarkerDict[egID] = markerKey

    # get all mouse markers
    results = db.sql('''select distinct a.accID as mgiId, m._Marker_key
	    from ACC_Accession a, MRK_Marker m
	    where a._MGIType_key = 2
	    and a._LogicalDB_key = 1
	    and a.prefixPart = 'MGI:'
	    and a._Object_key = m._Marker_key
	    and m._Marker_Status_key = 1
	    and m._Organism_key = 1''', 'auto')

    # removed per Richard
    # and a.preferred = 1
    for r in results:
	    mgiId = r['mgiId']
	    markerKey = r['_Marker_key']
	    mgiToMarkerDict[mgiId] = markerKey

    return

def parseFile():
    # Purpose: parse file into dictionary because new file has one line/MGI ID
    # Returns: 0
    # Assumes: humanEgToMouseMgiDict has been initialized
    # Effects: Reads file in file system
    # Throws: Nothing

    global humanEgToMouseMgiDict, hgncIdOnlyCount, mgiIdOnlyCount, egIdOnlyCount

    # These will be ignored by the load
    hgncIdOnlyCount = 0

    # These will be ignored by the load until we determine how to
    # cluster them in the absence of mouse MGI IDs
    mgiIdOnlyCount = 0

    # These will be ignored by the load
    egIdOnlyCount = 0

    # ignore header
    header = fpInFile.readline()

    for line in fpInFile.readlines():
	# 6/30 - file is no longer one line per cluster
	(hgncID, mgiID, egID) =  map(string.strip, string.split(line, TAB))
	#print 'egID: %s mgiID: %s' % (egID, mgiID) 

	# if both egID and mgiId columns are blank, skip and don't
        # report
        if egID == '' and mgiID == '':
            hgncIdOnlyCount += 1
            continue

        # if mgiId, and no egID skip and don't report
        elif egID == '':
            mgiIdOnlyCount += 1
            continue
	# if no mouse homology add to the count and mgiID value should be 'None'
        if mgiID == '':
             mgiID = 'None'
             egIdOnlyCount += 1
	# value  of humanEgToMouseMgiDict will be empty list if no mouse homology
	if not egID in humanEgToMouseMgiDict:
	    humanEgToMouseMgiDict[egID] = []
	# add the mouse homology to the dictionary
	humanEgToMouseMgiDict[egID].append(mgiID)
	
    return

def process():
    # Purpose: Create load ready file from HGNC file and the database
    # Returns: 0
    # Assumes: All lookup structures have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    global rptOne, rptTwo, rptDebug

    # parse the file into a data structure
    parseFile()

    # dictionary of id pairs to send to the clusterizer
    toClusterList = []
    lineCt = 0

    # for reporting - the actual line in the file including the header
    lineCt =  1  

    for egID in humanEgToMouseMgiDict:
	mgiIDList = humanEgToMouseMgiDict[egID]
	lineCt += 1
	
	# 1 means error on this line
	error = 0

	clusterFileLine = ''

	# current cluster - if there are no errors it will be added to
        # 'toClusterList'
        currentClusterList = []
	# report and skip lines where egID not in the database
	if egID and egID not in egToMarkerDict.keys():
	    toReport = '%s%s%s%s' % (egID, TAB, string.join(mgiIDList), CRT)
	    rptOne = '%s%s%s%s' % (rptOne, lineCt, TAB, toReport)
	    # if egID not in database continue to next input line
	    continue

	#
	# egID is in the database; check the mgi IDs
	#

	# add clusters with mouse to the list
	for id in mgiIDList:
	    id = string.strip(id)
	    # report and skip lines with mgiId not in the database
	    if id != 'None' and id not in mgiToMarkerDict.keys(): 
		error = 1
	        toReport = '%s%s%s%s' % (egID, TAB, id, CRT)
		rptTwo = '%s%s%s%s' % (rptTwo, lineCt, TAB, toReport)
		# No need to check any more ids, get out of the loop
		break
	    else:
		currentClusterList.append([egID, id])
	    clusterFileLine = ('%s%s%s%s%s' % \
		(clusterFileLine, egID, TAB, id, CRT))

	# if any mgi IDs not in database continue to next input line
	if error == 1:
		continue

	# if we have a cluster add it to the cluster list and to the file
	if currentClusterList != []:
	    # no errors so append the next cluster
	    toClusterList = toClusterList + currentClusterList

	    fpClustererFile.write(clusterFileLine)
	    # if we get here, we the egID is in the database and ALL the
	    # mgiIds are in the database
		
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
	# write debug to qc rpt
	rptDebug = '%s%s%s%s%s%s%s' % (rptDebug, idTuple, TAB, humanKeyList, TAB, mouseKeyList, CRT)
	fpLoadFile.write('%s%s%s%s' % (clusterId, TAB, keyString, CRT))

    return

def writeReports():
    # Purpose: writes out all sections of the QC report
    # Returns: 0
    # Assumes: rptOneand rptTwo have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    fpQcRpt.write('Count of lines with only an HGNC ID: %s (skipped)' % hgncIdOnlyCount)
    fpQcRpt.write('\nCount of lines with MGI ID, but no EG ID: %s (skipped)' % mgiIdOnlyCount)
    fpQcRpt.write('\nCount of lines with EG ID, but no MGI ID: %s (loaded; human singletons)\n\n\n' % egIdOnlyCount)

    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)

    fpQcRpt.write(rptDebug)
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
