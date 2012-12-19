#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       From HomoloGene input file create load ready file
#
# Usage: createInputFiles.py
#
# Inputs:
#	1. HomoloGene file tab-delimited in following format:
#	    1. 
#	    2. 
#	    3. 
#           4.
#           5.
#           5.
#	2. Configuration - see homologyload.config
#	  1.
#         2.
#         3.
#         4.
#         5.
# Outputs:
#	 1. tab delimited file:
#           1. 
#           2. 
#           3. 
#           4. 
#           5. 
#           6. 
#           7. 
#           8. 
#	 2. log file
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


# constants
TAB= '\t'
CRT = '\n'
SPACE = ' '

ids = os.environ['LOADED_TAXIDS']
taxIds = ids.split(',')
taxIdList = []

# remove any spaces
for id in taxIds:
    taxIdList.append(id.strip())

# From the database: egID:[list of Marker instances], ...}
egToMarkerDict = {}

# From the database: _Marker_key:[list of egIds], ...}
mrkToMultiEGDict = {}

# egIDs mapped to their inputlines
# {egId:[list of input lines], ...}
egIdToLineDict = {}

# paths to input and output files
inFilePath = os.environ['INPUT_FILE']

# This is the cleaned up input file
outFilePath = os.environ['INPUT_FILE_LOAD']

# The QC report path
qcRptPath = os.environ['QC_RPT']

sep = '--------------------------------------------------\n'
head = 'hgID%staxID%segID%s' % (TAB, TAB, CRT)
rptOne = 'EG IDs found in > 1 class in the input file\n\n%s%s' % (head, sep)

rptTwo = '\n\nEG IDs in input file and not in MGI, therefore will not participate in a class\n\n%s%s' % (head, sep)

head2 = 'hgID%staxID%segID%sOther Markers%s' % (TAB, TAB, TAB, CRT)
rptThree = '\n\nEG IDs in the input file associated with >1 marker in MGI\n\n%s%s' % (head2, sep)

head3 = 'hgID%staxID%segID%sOther EG IDs%s' % (TAB, TAB, TAB, CRT)
rptFour = '\n\nMarkers associated with input EG ID and also associated with other EG IDs in MGI\n\n%s%s' % (head3, sep)



# file descriptors

fpInFile = ''
fpOutFile = ''
fpSanityRpt = ''
fpQcRpt = ''

class Marker:
    def __init__(self, markerKey, egId, symbol, organism):
	self.k = markerKey
	self.e = egId
	self.s = symbol
	self.o = organism

user = os.environ['MGD_DBUSER']
passwordFileName = os.environ['MGD_DBPASSWORDFILE']
db.useOneConnection(1)
db.set_sqlUser(user)
db.set_sqlPasswordFromFile(passwordFileName)

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
#
# Initialize
#

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

#
# Process
#

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

print 'iterating over input file lookup'
for egId in egIdToLineDict.keys():
    lineList = egIdToLineDict[egId]

    # egId participates in > class in the input file
    if len(lineList) > 1:
	rptOne = rptOne + 'egID: %s\n' % (egId)
	for line in lineList:
	    rptOne = rptOne + '    %s\n' % \
		(string.strip(string.join(line, ',')).split('\t'))

    # we have only one line, format it
    line = (string.strip(string.join(lineList, ',')).split('\t'))

    # egId not associated with a marker in MGI
    if not egToMarkerDict.has_key(egId):
	rptTwo = rptTwo + '%s%s%s%s%s%s' % \
	    (line[0], TAB, line[1], TAB, egId, CRT)

    else:
	markerList = egToMarkerDict[egId]
	if len(markerList) > 1:	# egId associated with > 1 marker in MGI
	    symbols = ''
	    for m in markerList:
		symbols = symbols + '%s|%s ' % (m.s, m.o)
	    rptThree = rptThree + '%s%s%s%s%s%s%s%s' % \
		(line[0], TAB, line[1], TAB, egId, TAB, symbols, CRT)
	else: # there is only one marker in the db for this egId
	    m = egToMarkerDict[egId][0]
	    # this marker is also associated with other egId(s)
	    if mrkToMultiEGDict.has_key(m.k):
		mDbList = mrkToMultiEGDict[m.k]
		otherEgIds = ''
		for mDb in mDbList:
                        otherEgIds = otherEgIds +  '%s ' % mDb.e
		rptFour = rptFour + '%s%s%s%s%s%s%s%s' % \
		    (line[0], TAB, line[1], TAB, egId, TAB, otherEgIds, CRT)
	    # write it out
	    else:
		#print 'line: %s' % line
		#print 'k:%s e:%s s:%s o:%s' % (m.k, m.e, m.s, m.o)
		fpOutFile.write('%s%s%s%s%s%s%s%s%s%s%s%s' % \
		    (line[0], TAB, line[1], TAB, egId, TAB, m.k, TAB, m.s, TAB, m.o, CRT))
    
# (hgId, taxId, egId, junk1, junk2, junk3)
# postprocess

fpQcRpt.write(rptOne)
fpQcRpt.write(rptTwo)
fpQcRpt.write(rptThree)
fpQcRpt.write(rptFour)

fpInFile.close()
fpOutFile.close()
fpQcRpt.close()
db.useOneConnection(0)
print '%s' % mgi_utils.date()
