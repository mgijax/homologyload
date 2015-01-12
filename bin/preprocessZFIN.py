#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       From ZFIN input files create load ready file
#
# Usage: preprocessZFIN.py
#
# Inputs:
#	1. ZFIN mouse orthology file tab delimited:
#	    1. ZFIN ID 	
# 	    2. ZFIN Symbol (not used)
# 	    3. ZFIN Name (not used)
# 	    4. Mouse Symbol (not used)
# 	    5. Mouse Name (not used)
# 	    6. MGI ID 
#	    7. Gene ID (not used)
#
#       2. ZFIN gene file tab delimited:
#           1. ZFIN ID
#           2. SO ID (not used)
#           3. Symbol (not used)
#           4. NCBI Gene ID
#
#	3. ZFIN Genes with Expression Assay Records file tab delimited:
# 	    1. ZFIN ID 	
# 	    2. Gene Symbol (not used)
# 	    3. EST ID (Optional) (not used)
# 	    4. EST Symbol (Optional) (not used)
# 	    5. Expression Type (not used)	
# 	    6. Expression ID (not used)
# 	    7. Publication ID (not used)
# 	    8. Genotype ID (not used)
# 	    9. Environment ID (not used)
# 	    10 Probe Quality (optiona 0 - 5 rating) (not used)
#
#	4. Configuration - see homologyload.config
#	  . INPUT_FILE_* - path to file in load input file directory
#         . INPUT_FILE_LOAD - path to load-ready file
#         . QC_RPT - path to QC report
#         . MGD_DBUSER - database user
#	  . MGD_DBPASSWORDFILE - database password
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
###########################################################################

import sys
import os
import string
import Set

import mgi_utils
import loadlib
import db
import clusterize

print '%s' % mgi_utils.date()


####################################
#
# Globals
#
####################################

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

# input files from ZFIN
inFileGenePath = os.environ['INPUT_FILE_GENE']
inFileOrthoPath = os.environ['INPUT_FILE_ORTHO']
inFileExprPath = os.environ['INPUT_FILE_EXPR']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

# ZFIN gene IDs from the expression file
exprSet = set([])

# ZFIN gene IDs mapped to their NCBI IDs from gene file
geneDict = {}

# XFIN gene IDs mapped to one or more mouse MGI IDs from orthology file
mouseDict = {}

#
# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'

rptOne = 'ZFIN IDs from expression file not in either gene or orthos file%s%s%s' % (CRT, CRT, sep)
rptTwo = '%s%sZFIN NCBI ID not in database%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptThree = '%s%sMouse MGI ID not in database%s%s%s' % (CRT, CRT,CRT, CRT, sep)

#
# file descriptors
#

fpGeneFile = ''
fpOrthoFile = ''
fpExprFile = ''

# idList replaces this file, but keeping for time being
fpLoadFile = ''
fpQcRpt = ''

#
# class to hold information about an MGI marker
#
class Marker:
    def __init__(self, markerKey, egId, symbol):
	self.k = markerKey
	self.e = egId
	self.s = symbol

def init():
    global egToMarkerDict, mgiToMarkerDict
    global fpGeneFile, fpOrthoFile, fpExprFile
    global fpLoadFile, fpQcRpt

    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    try:
	fpGeneFile = open(inFileGenePath, 'r')
    except:
	exit('Could not open file for reading %s\n' % inFileGenePath)

    try:
        fpOrthoFile = open(inFileOrthoPath, 'r')
    except:
        exit('Could not open file for reading %s\n' % inFileOrthoPath)

    try:
        fpExprFile = open(inFileExprPath, 'r')
    except:
        exit('Could not open file for reading %s\n' % inFileExprPath)

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
	and m._Organism_key = 84''', 'auto')
    #
    # create Human egID to marker lookup from database
    #
    #print 'creating egID to marker lookup'
    for r in results:
	egId = r['egId']
	markerKey = r['_Marker_key']
	#print 'egId: %s' % egId
	egToMarkerDict[egId] = markerKey
    #print ' size of egToMarkerDict %s' % len(egToMarkerDict)

    # get all mouse markers
    #print 'creating mouse marker lookup'
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
    #print 'size of mgiToMarkerDict %s' % len(mgiToMarkerDict)

    return

def processInputFiles():
    global exprSet, geneDict, mouseDict
    for line in fpExprFile.readlines():
	tokens = string.split(line, TAB)
	zfinID = string.strip(tokens[0])
	if zfinID.startswith('ZDB-GENE'):
	    exprSet.add(zfinID)
    #print exprSet
    for line in fpGeneFile.readlines():
	tokens = string.split(line, TAB)
	zfinID = string.strip(tokens[0])
	if zfinID.startswith('ZDB-GENE'):
	    ncbiID = string.strip(tokens[3])
	    geneDict[zfinID] = ncbiID
    #print geneDict.keys()
    for line in fpOrthoFile.readlines():
	tokens = string.split(line, TAB)
	zfinID = string.strip(tokens[0])
	if zfinID.startswith('ZDB-GENE'):
	    mgiID = string.strip(tokens[5])
	    if not mouseDict.has_key(zfinID):
		mouseDict[zfinID] = []
	    mouseDict[zfinID].append(mgiID)
    #print mouseDict.keys()

def process():
    global rptOne, rptTwo, rptThree
    # egToMarkerDict, mgiToMarkerDict

    # dictionary of id pairs to send to the clusterizer
    toClusterList = []
    print 'processing input file'
    zfinIdNotInSet = set([])
    ncbiNotInDBSet = set([])
    mgiNotInDBSet = set([])

    for zfinID in exprSet:
	# Join to gene and orthos to get ncbi and mgi IDs
	print zfinID
	ncbiID = ''
	mgiID = ''
	error = 0
	if zfinID in geneDict.keys():
	    ncbiID = geneDict[zfinID]
	else:
	    error = 1
	    zfinIdNotInSet.add(zfinID)
	if zfinID in mouseDict.keys():
	    mgiIdList = mouseDict[zfinID]
	else:
	    error = 1
	    zfinIdNotInSet.add(zfinID)
	if error:
	    print 'continuing: %s' % zfinID
	    continue
	# verify zfin NCBI ID in database
	if ncbiID not in egToMarkerDict.keys():
	    ncbiNotInDBSet.add(ncbiID)
	    error = 1
	# verify mouse mgiIDs in database
	currentClusterList = []
	for mgiID in mgiIdList:
	    currentClusterList.append([ncbiID, mgiID])
	    if mgiID not in mgiToMarkerDict.keys():
		mgiNotInDBSet.add(mgiID)
		error = 1
	if error:
	    print 'ncbi or mgi not in db continuing: %s' % zfinID
	    continue
	else:
	    # no errors so append the next cluster
	    toClusterList = toClusterList + currentClusterList
    clusterDict = clusterize.cluster(toClusterList, 'ZFIN')
    # now resolve the ids to database keys; zfin and mouse gene keys
    for clusterId in clusterDict.keys():
        idTuple = clusterDict[clusterId]
        zfinKeyList = []
        mouseKeyList = []
        for id in idTuple:
            if id.startswith('MGI:'):
                mouseKeyList.append(str(mgiToMarkerDict[id]))
            else:
                zfinKeyList.append(str(egToMarkerDict[id]))
        # we want mouse to come before zfin for cluster member sequence numbering
        keyList = mouseKeyList + zfinKeyList
        keyString = ', '.join(keyList)
        fpLoadFile.write('%s%s%s%s' % (clusterId, TAB, keyString, CRT))


    for id in zfinIdNotInSet:
	rptOne = '%s%s%s' % (rptOne, id, CRT)
    rptOne = '%s%sTotal IDs: %s%s' % (rptOne, CRT, len(zfinIdNotInSet), CRT) 

    for id in ncbiNotInDBSet:
	rptTwo =  '%s%s%s' % (rptTwo, id, CRT)
    rptTwo = '%s%sTotal IDs: %s%s' % (rptTwo, CRT, len(ncbiNotInDBSet), CRT)

    for id in mgiNotInDBSet:
	rptThree = '%s%s%s' % (rptThree, id, CRT)
    rptThree = '%s%sTotal IDs: %s%s' % (rptThree, CRT, len(mgiNotInDBSet), CRT)

    #print 'zfinID: %s, ncbiID: %s, mgiID: %s' % (zfinID, ncbiID, mgiID)
    return

def writeReports():
    #print 'writing reports'
    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)
    fpQcRpt.write(rptThree)

    return

def closeFiles():
    #print 'closing files'
    fpGeneFile.close()
    fpOrthoFile.close()
    fpExprFile.close()
    fpLoadFile.close()
    fpQcRpt.close()

    # close the database connection
    db.useOneConnection(0)
    
    return
init()
processInputFiles()
process()
writeReports()
closeFiles()

print '%s' % mgi_utils.date()
