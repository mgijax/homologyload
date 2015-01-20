#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       From GEISHA (chicken database) input files create load ready file
#
# Usage: preprocessGEISHA.py
#
# Inputs:
#	1. GEISHA mouse orthology file tab delimited:
#	    1. Chicken NCBI ID
# 	    2. Chicken gene name (not used)
#	    3. GEISHA ID  (not used)
#	    4. Human NCBI gene ID       
# 	    5. Mouse NCBI ID
#	    6. Xenopus NCBI gene ID (not used)
#	    7. Zebrafish NCBI Gene ID (not used)
#
#	2. GEISHA Genes with Expression Assay Records file tab delimited:
#	    1. Chicken NCBI gene ID     
# 	    2. Chicken gene name        
#	    3. GEISHA ID BirdBase ID    
#	    4. Ensembl ID       
#	    5. GO ids   
#	    6. stages   
#	    7. locations
#
#	3. Configuration - see homologyload.config
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

# EG ID/Chicken Marker associations from the database
# {ncbiID:chicken marker key, ...}
ncbiToChickenDict = {}

# EG ID/Mouse Marker associations from the database
# {ncbiID:mouse marker key, ...}
ncbiToMouseDict = {}

#
# paths to input and output files
#

# orthology input file
inFileOrthoPath = os.environ['INPUT_FILE_ORTHO']

# The expression input file
inFileExprPath = os.environ['INPUT_FILE']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

# Chicken NCBI gene IDs from the expression file
exprSet = set([])

# Chicken NCBI gene IDs mapped to one or more mouse NCBI IDs from orthology file
mouseDict = {}

#
# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'

rptOne = 'Chicken NCBI IDs from expression file not in orthos file%s%s%s' % (CRT, CRT, sep)
rptTwo = '%s%sChicken NCBI IDs not in database%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptThree = '%s%sMouse NCBI IDs not in database%s%s%s' % (CRT, CRT,CRT, CRT, sep)

#
# file descriptors
#

fpOrthoFile = ''
fpExprFile = ''

# idList replaces this file, but keeping for time being
fpLoadFile = ''
fpQcRpt = ''

def init():
    global ncbiToChickenDict, ncbiToMouseDict
    global fpOrthoFile, fpExprFile
    global fpLoadFile, fpQcRpt

    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

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


    # get all chicken markers that are associated with ncbi IDs
    results = db.sql('''select distinct a.accid as ncbiId, m._Marker_key
	from ACC_Accession a, MRK_Marker m
	where a._MGIType_key = 2
	and a._LogicalDB_key = 55
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Status_key in (1,3)
	and m._Organism_key = 63''', 'auto')
    #
    # create Chicken ncbiID to marker lookup from database
    #
    #print 'creating chicken ncbiID to marker lookup'
    for r in results:
	ncbiId = r['ncbiId']
	markerKey = r['_Marker_key']
	#print 'ncbiId: %s' % ncbiId
	ncbiToChickenDict[ncbiId] = markerKey
    #print ' size of ncbiToChickenDict %s' % len(ncbiToChickenDict)

    # get all mouse markers
    #print 'creating mouse marker lookup'
    results = db.sql('''select distinct a.accID as ncbiId, m._Marker_key
	    from ACC_Accession a, MRK_Marker m
	    where a._MGIType_key = 2
	    and a._LogicalDB_key = 55
	    and a._Object_key = m._Marker_key
	    and m._Marker_Status_key in (1,3)
	    and m._Organism_key = 1''', 'auto')

    for r in results:
	    ncbiId = r['ncbiId']
	    markerKey = r['_Marker_key']
	    ncbiToMouseDict[ncbiId] = markerKey
    #print 'size of ncbiToMouseDict %s' % len(ncbiToMouseDict)

    return

def processInputFiles():
    global exprSet, mouseDict

    # remove header
    header = fpExprFile.readline()
    line =  fpExprFile.readline()
    while line:
	tokens = string.split(line, TAB)
	ncbiChickenID = string.strip(tokens[0])
	exprSet.add(ncbiChickenID)
	line =  fpExprFile.readline()
    #print exprSet
    # remove header
    header = fpOrthoFile.readline()
    line = fpOrthoFile.readline()
    while line:
	tokens = string.split(line, TAB)
	ncbiChickenID = string.strip(tokens[0])
	ncbiMouseIDs = string.strip(tokens[4])
	# skip if no mouse ID(s)  on this line
	if ncbiMouseIDs == '':
	    line = fpOrthoFile.readline()
	    continue
	else:
	    print ncbiMouseIDs
	    print string.split(ncbiMouseIDs, ',')
	    mouseList = string.split(ncbiMouseIDs, ',')
	    mouseDict[ncbiChickenID] = mouseList
	#for id in mouseList:
	#    if not mouseDict.has_key(ncbiChickenID):
	#	mouseDict[ncbiChickenID] = []
	#    mouseDict[ncbiChickenID].append(id)
	line = fpOrthoFile.readline()
    #print mouseDict.keys()

def process():
    global rptOne, rptTwo, rptThree
    # dictionary of id pairs to send to the clusterizer
    toClusterList = []
    print 'processing input file'
    chickenIdNotInSet = set([])
    chickenIdNotInDBSet = set([])
    mouseIdNotInDBSet = set([])

    for chickenID in exprSet:
	# Join to orthos to get mouse ncbi Ids
	#print ncbiChickenID
	mouseID = ''
	error = 0
	# get the mouse orthology, if none report
	if chickenID in mouseDict.keys():
	    mouseIdList = mouseDict[chickenID]
	else:
	    error = 1
	    chickenIdNotInSet.add(chickenID)
	if error:
	    continue
	# verify chicken NCBI ID in database
	if chickenID not in ncbiToChickenDict.keys():
	    chickenIdNotInDBSet.add(chickenID)
	    error = 1
	# verify mouse NCBI ID in database
	currentClusterList = []
	for mouseID in mouseIdList:
	    currentClusterList.append([chickenID, mouseID])
	    if mouseID not in ncbiToMouseDict.keys():
		mouseIdNotInDBSet.add(mouseID)
		error = 1
	if error:
	    continue
	else:
	    # no errors so append the next cluster
	    toClusterList = toClusterList + currentClusterList
    clusterDict = clusterize.cluster(toClusterList, 'GEISHA')
    # now resolve the ids to database keys; chicken and mouse gene keys
    # NOT SORTING BY ORGANISM for sequenceNum because we don't need to. 
    # If we find we need for this load we will need to create a mouse and a 
    # chicken lookup by NCBI ID to determine which organism in order to 
    # order correctly (mouse first then chicken)
    for clusterId in clusterDict.keys():
        idTuple = clusterDict[clusterId]
	mouseKeyList = []
	chickenKeyList = []
	for id in idTuple:
	    if id in ncbiToMouseDict.keys():
		mouseKeyList.append(str(ncbiToMouseDict[id]))
	    elif id in ncbiToChickenDict.keys():
		chickenKeyList.append(str(ncbiToChickenDict[id]))
	    else:
		print 'not chicken or mouse'
	keyList = mouseKeyList + chickenKeyList
	keyString = ', '.join(keyList)
        fpLoadFile.write('%s%s%s%s' % (clusterId, TAB, keyString, CRT))

    for id in chickenIdNotInSet:
	rptOne = '%s%s%s' % (rptOne, id, CRT)
    rptOne = '%s%sTotal IDs: %s%s' % (rptOne, CRT, len(chickenIdNotInSet), CRT) 

    for id in chickenIdNotInDBSet:
	rptTwo =  '%s%s%s' % (rptTwo, id, CRT)
    rptTwo = '%s%sTotal IDs: %s%s' % (rptTwo, CRT, len(chickenIdNotInDBSet), CRT)

    for id in mouseIdNotInDBSet:
        rptThree =  '%s%s%s' % (rptThree, id, CRT)
    rptThree = '%s%sTotal IDs: %s%s' % (rptThree, CRT, len(mouseIdNotInDBSet), CRT)

    return

def writeReports():
    #print 'writing reports'
    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)
    fpQcRpt.write(rptThree)

    return

def closeFiles():
    #print 'closing files'
    fpOrthoFile.close()
    fpExprFile.close()
    fpLoadFile.close()
    fpQcRpt.close()

    # close the database connection
    db.useOneConnection(0)
    
    return

print 'initializing'
init()
print 'processInputFiles'
processInputFiles()
print 'processing clusters'
process()
print 'writing reports'
writeReports()
print 'closing files'
closeFiles()

print '%s' % mgi_utils.date()
