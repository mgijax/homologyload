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
#	2. Configuration - see homologyload.config
#	  1. 
#         2. INPUT_FILE - path to input file
#         3. INPUT_FILE_LOAD - path to cleaned up load-ready file
#         4. QC_RPT - path to QC report
#         5. MGD_DBUSER - database user
#	  6. MGD_DBPASSWORDFILE - database password
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
# {egID:[list of Marker instances], ...}
egToMarkerDict = {}

# list of all preferred mouse MGI IDs in the database
mouseMarkerIdList = []
#
# paths to input and output files
#

# input file from Homologene
inFilePath = os.environ['INPUT_FILE']

# This is the cleaned up load-ready input file
clustererFilePath = os.environ['INPUT_FILE_CLUSTERER']

# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'
head = 'Line#%sHuman egID%sMouse MGI IDs%s HGNC ID%s' % (TAB, TAB, TAB, CRT)

rptOne = 'Lines where Human egID not in database%s%s%s%s' % (CRT, CRT, head, sep)

rptTwo = '%s%sLines where a Mouse MGI ID not in database%s%s%s%s' % (CRT, CRT,CRT, CRT, head, sep)

#
# file descriptors
#

fpInFile = ''
fpClustererFile = ''
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
    global egToMarkerDict, mouseMarkerIdList
    global fpInFile, fpClustererFile, fpQcRpt

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
	fpQcRpt = open(qcRptPath, 'w')
    except:
	exit('Could not open file for writing %s\n' % qcRptPath)


    # get all human markers that are associated with egIds
    results = db.sql('''select distinct a.accid as egId, m._Marker_key, m.symbol
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
    print 'creating egID to marker lookup'
    for r in results:
	egId = r['egId']
	markerKey = r['_Marker_key']
	symbol = r['symbol']
	#print 'egId: %s' % egId
	if not egToMarkerDict.has_key(egId):
	    egToMarkerDict[egId] = []
	    egToMarkerDict[egId].append( Marker(markerKey, egId, symbol) )
    print ' size of egToMarkerDict %s' % len(egToMarkerDict)

    # get all mouse markers
    print 'creating mouse marker lookup'
    results = db.sql('''select distinct a.accID as mgiID
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
	    mouseMarkerIdList.append(r['mgiID'])
    print 'size of mouseMarkerIdList %s' % len(mouseMarkerIdList)

    return

def process():
    global rptOne, rptTwo, hgncIdOnlyCount, mgiIdOnlyCount

    print 'processing input file'
    lineCt = 0
    # These will be ignored by the load
    hgncIdOnlyCount = 0

    # These will be ignored by the load until we determine how to 
    # cluster them in the absence of mouse MGI IDs
    mgiIdOnlyCount = 0

    for line in fpInFile.readlines():
	#print 'line: %s' % line
	lineCt += 1
	(egId, mgiIDs, junk1) = string.split(line[:-1], TAB)
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
	#print 'mgiIDs: %s' % mgiIDs
	egId = string.strip(egId)
	error = 0
	clusterFileLine = ''
	if egId and egId not in egToMarkerDict.keys():
	    error = 1
	    rptOne = '%s%s%s%s' % (rptOne, lineCt, TAB, line)
	    print 'writing to rptOne: %s %s' % (lineCt, line)
	mgiIDList = map(string.strip, string.split(mgiIDs, ','))
	# line to write to the cluster
	for id in mgiIDList:
	    clusterFileLine = ('%s%s%s%s%s' % \
		(clusterFileLine, egId, TAB, id, CRT))
	    if id != 'None' and id not in mouseMarkerIdList:
		error = 1
		print 'id not in MGI: %s' % id
		rptTwo = '%s%s%s%s' % (rptTwo, lineCt, TAB, line)
		print 'writing to rptTwo: %s %s' % (lineCt, line)
		break
	if not error:
	    fpClustererFile.write(clusterFileLine)
    fpClustererFile.close()
    clusterDict = clusterize.cluster(clustererFilePath, 'HGNC')
    print 'clusterDict: %s' % clusterDict
    return

def writeReports():
    print 'writing reports'
    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)

    return

def cluster():
    fpCF = open(clustererFilePath, 'r')
    humanToMouseDict = {}
    mouseToHumanDict = {}
    for line in fpCF.readlines():
	humanId, mouseId = string.split(line)
	#print '%s, %s' % (humanId, mouseId)
	if not humanToMouseDict.has_key(humanId):
	    humanToMouseDict[humanId] = []
	humanToMouseDict[humanId].append(mouseId)
	if mouseId != 'None':
	    if not mouseToHumanDict.has_key(mouseId):
		mouseToHumanDict[mouseId] = []
	    mouseToHumanDict[mouseId].append(humanId)
    print 'humanToMouseDict'
    print humanToMouseDict
    print 'mouseToHumanDict'
    print mouseToHumanDict
    fpCF.close()

    # set of all clusters, uniq set of tuples
    clusterSet = set()
    clusterCt = 0 
    for humanId in humanToMouseDict.keys():
	print 'One: ' + humanId
	# current cluster
	currentList = []
	#currentList.append(humanId)
	for mouseId in humanToMouseDict[humanId]:
	    print 'Two: ' + mouseId
	    if mouseId != 'None':
		currentList.append(mouseId)
	    if mouseToHumanDict.has_key(mouseId):
		for hId in mouseToHumanDict[mouseId]:
		    print 'Three: ' + hId
		    if hId not in currentList:
			currentList.append(hId)
	print 'Four: %s' % currentList
	currentList.sort()
	clusterSet.add(tuple(currentList))
    fpLF = open(os.environ['INPUT_FILE_LOAD'], 'w')
    clusterCt = 0
    for c in clusterSet:
  	clusterCt +=1
	nextId = 'HGNC:%s' % clusterCt
	fpLF.write('%s%s%s%s' % (nextId, TAB, ', '.join(c), CRT))
    fpLF.close()


def closeFiles():
    print 'closing files'
    fpInFile.close()
    fpQcRpt.close()

    # close the database connection
    db.useOneConnection(0)
    
    return
init()
process()
writeReports()
closeFiles()

print 'lines with HGNC ID only: %s' % hgncIdOnlyCount
print 'lines no egID: %s' % mgiIdOnlyCount

print '%s' % mgi_utils.date()
