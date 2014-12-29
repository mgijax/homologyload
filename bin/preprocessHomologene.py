#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       From HomoloGene input file create load ready file
#
# Usage: runQC.py
#
# Inputs:
#	1. HomoloGene file tab-delimited in following format:
#	    1. Cluster ID - identifies the Homologene cluster
#	    2. Taxonomy ID - identifies the organism for this marker
#	    3. EntrezGene ID - identifies the marker itself
#           4. Marker symbol
#           5. MGI ID - NCBI identifier for the RefSeq
#           6. Protein sequence ID - RefSeq
#
#	2. Configuration - see homologyload.config
#	  1. LOADED_TAXIDS - list of organisms we are loadin
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

import mgi_utils
import loadlib
import db

print '%s' % mgi_utils.date()


####################################
#
# Globals
#
####################################

# constants
TAB= '\t'
CRT = '\n'

#
# get the tax IDs we are loading homologies for
#
ids = os.environ['LOADED_TAXIDS']
taxIds = ids.split(',')
taxIdList = []
for id in taxIds:	# remove any spaces
    taxIdList.append(id.strip())


# EG ID/Marker associations from the database
# {egID:[list of Marker instances], ...}
egToMarkerDict = {}

# Marker/EG ID associations from the database
# { _Marker_key:[list of Marker instances], ...}
mrkToMultiEGDict = {}

# egIDs from the input mapped to their input lines
# {egId:[list of input lines], ...}
egIdToLineDict = {}

#
# paths to input and output files
#

# input file from Homologene
inFilePath = os.environ['INPUT_FILE']

# This is the cleaned up load-ready input file
outFilePath = os.environ['INPUT_FILE_LOAD']

# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'
head1 = 'hgID%staxID%segID%s' % (TAB, TAB, CRT)
head2 = 'hgID%staxID%segID%sOther Markers%s' % (TAB, TAB, TAB, CRT)
head3 = 'hgID%staxID%segID%sOther EG IDs%s' % (TAB, TAB, TAB, CRT)

rptOne = 'EG IDs found in > 1 class in the input file\n\n%s%s' % (head1, sep)

rptTwo = '\n\nEG IDs in input file and not in MGI, therefore will not participate in a class\n\n%s%s' % (head1, sep)

rptThree = '\n\nEG IDs in the input file associated with >1 marker in MGI\n\n%s%s' % (head2, sep)

rptFour = '\n\nMarkers associated with input EG ID and also associated with other EG IDs in MGI\n\n%s%s' % (head3, sep)

#
# file descriptors
#

fpInFile = ''
fpOutFile = ''
fpQcRpt = ''

#
# class to hold information about an MGI marker
#
class Marker:
    def __init__(self, markerKey, egId, symbol, organism):
	self.k = markerKey
	self.e = egId
	self.s = symbol
	self.o = organism

####################################
#
# create database connection
#
####################################

user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
db.useOneConnection(1)
db.set_sqlUser(user)
db.set_sqlPasswordFromFile(passwordFileName)

####################################
#
# Initialize file descriptors
#
####################################

try:
    fpInFile = open(inFilePath, 'r')
except:
    exit('Could not open file for reading %s\n' % inFilePath)

try:
    fpOutFile = open(outFilePath, 'w')
except:
    exit('Could not open file for writing %s\n' % outFilePath)

try:
    fpQcRpt = open(qcRptPath, 'w')
except:
    exit('Could not open file for writing %s\n' % qcRptPath)


####################################
#
# Load all Lookup dictionaries
#
####################################

#
# create mrkToMultiEGDict
#

# get all markers that are associated with egIds, all organisms
db.sql('''select distinct a.accid as egId, m._Marker_key, m.symbol, o.commonName
	into #eg
	from ACC_Accession a, MRK_Marker m, MGI_Organism o
	where a._MGIType_key = 2
	and a._LogicalDB_key = 55
	and a.preferred = 1
	and a._Object_key = m._Marker_key
	and m._Marker_Status_key in (1,3)
	and m._Organism_key = o._Organism_key''', None)

results = db.sql('''select * from #eg
	group by _Marker_key having count(*) > 1''', 'auto')

print 'creating marker to multi egID lookup'
for r in results:
    egId = r['egId']
    markerKey = r['_Marker_key']
    symbol = r['symbol']
    organism = r['commonName']
    if not mrkToMultiEGDict.has_key(markerKey):
	mrkToMultiEGDict[markerKey] = []
    mrkToMultiEGDict[markerKey].append( Marker(markerKey, egId, symbol, organism) )

#
# create egToMarkerDict
#

results = db.sql('''select * from #eg''', 'auto')

print 'creating egID to marker lookup'
for r in results:
    egId = r['egId']
    markerKey = r['_Marker_key']
    symbol = r['symbol']
    organism = r['commonName']
    #print 'egId: %s' % egId
    if not egToMarkerDict.has_key(egId):
	egToMarkerDict[egId] = []
    egToMarkerDict[egId].append( Marker(markerKey, egId, symbol, organism) )

print 'creating input file lookup'
for line in fpInFile.readlines():

    (hgId, taxId, egId, junk1, junk2, junk3) = string.split(line[:-1], TAB)

    hgId = string.strip(hgId)
    taxId = string.strip(taxId)
    egId = string.strip(egId)
    if taxId not in taxIdList:
        continue
    if not egIdToLineDict.has_key(egId):
        egIdToLineDict[egId] = []
    egIdToLineDict[egId].append(line)

#
# Process
#

print 'iterating over input file lookup'
for egId in egIdToLineDict.keys():
    lineList = egIdToLineDict[egId]
    
    #
    # if EG IDs found in > 1 class in the input file - report and skip 
    # 
    if len(lineList) > 1:
	rptOne = rptOne + 'egID: %s\n' % (egId)
	for line in lineList:
	    rptOne = rptOne + '    %s\n' % \
		(string.strip(string.join(line, ',')).split('\t'))
	continue

    # we have only one class, format  the line
    line = (string.strip(string.join(lineList, ',')).split('\t'))

    #
    # EG ID not in MGI, therefore will not participate in a class
    # Report. All other members of class will be loaded
    # 
    if not egToMarkerDict.has_key(egId):
	rptTwo = rptTwo + '%s%s%s%s%s%s' % \
	    (line[0], TAB, line[1], TAB, egId, CRT)

    else: # EG ID in MGI and not found in > 1 class in the input
	markerList = egToMarkerDict[egId]
	#
	# EG ID in the input associated with > 1 marker in MGI - report and skip
	#
	if len(markerList) > 1:	# egId associated with > 1 marker in MGI
	    symbols = ''
	    for m in markerList:
		symbols = symbols + '%s|%s ' % (m.s, m.o)
	    rptThree = rptThree + '%s%s%s%s%s%s%s%s' % \
		(line[0], TAB, line[1], TAB, egId, TAB, symbols, CRT)
	else: # there is only one marker in the db for this egId
	    m = egToMarkerDict[egId][0]
	    #
	    # Markers associated with input EG ID and also associated with 
	    # other EG IDs in MGI - report and skip
	    #
	    if mrkToMultiEGDict.has_key(m.k):
		mDbList = mrkToMultiEGDict[m.k]
		otherEgIds = ''
		for mDb in mDbList:
                        otherEgIds = otherEgIds +  '%s ' % mDb.e
		rptFour = rptFour + '%s%s%s%s%s%s%s%s' % \
		    (line[0], TAB, line[1], TAB, egId, TAB, otherEgIds, CRT)
	    # All QC checks passed - write out to the load-ready file
	    else:
		#print 'line: %s' % line
		#print 'k:%s e:%s s:%s o:%s' % (m.k, m.e, m.s, m.o)
		hgId = line[0]
		taxId = line[1]
		symbol =  m.s
		markerKey = m.k
		organism = m.o
		fpOutFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' % \
		    ( organism, TAB, symbol, TAB, hgId, TAB, egId, TAB, markerKey, TAB, taxId, CRT))
   
#####################################
# 
# postprocess
#
#####################################

#
# write out all QC reports
#
fpQcRpt.write(rptOne)
fpQcRpt.write(rptTwo)
fpQcRpt.write(rptThree)
fpQcRpt.write(rptFour)

#
# close all file descriptors
#
fpInFile.close()
fpOutFile.close()
fpQcRpt.close()

# close the database connection
db.useOneConnection(0)

print '%s' % mgi_utils.date()