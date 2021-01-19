##########################################################################
#
# Purpose:
#       From HGNC input file create load ready file
#
# Usage: preprocessHGNC.py
#
# Inputs:
#       1. Alliance file tab-delimited in following format:
#
#           1. Gene1ID 
#           2. Gene1Symbol 
#           3. Gene1SpeciesTaxonID 
#           4. Gene1SpeciesName 
#           5. Gene2ID
#           6. Gene2Symbol
#           7. Gene2SpeciesTaxonID
#           8. Gene2SpeciesName
#           9-13 (not used)
#
#       2. Configuration - see alliance_directload.config
#
# Outputs:
#        1. load ready file
#        2. QC report file
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
# sc   01/15/2021
#       - initial implementation
###########################################################################

import copy
import os
import string
import mgi_utils
import clusterize
import db

###--- globals ---###

# constants
TAB= '\t'
CRT = '\n'

# [ 'human', 'mouse, laboratory', 'rat', 'zebrafish' ]

organismOrder = [2, 1, 40, 84]

#
# paths to input and output files
#

# input file from HGNC
inFilePath = os.environ['INPUT_FILE']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

# The QC report
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'

rptOne = 'Lines where a Mouse MGI ID not in database %s%s%s' % (CRT, sep, CRT)
rptOne = rptOne + 'LineNum%sline%s' % (TAB, CRT)
rptTwo = '%s%sLines where a Homology ID not in database%s%s%s%s' % (CRT, CRT, CRT, CRT, sep, CRT)
rptTwo = rptTwo + 'LineNum%sline%s' % (TAB, CRT)

#
# file descriptors
#

fpInFile = ''
fpLoadFile = ''
fpQcRpt = ''

#
# Database Lookups
#

# Rat, Human Zebra Fish lookup
# ( ID:_Marker_key)
homologyLookup = {}

# Mouse gene lookup
# (ID:_Marker_key)
mouseLookup = {}

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
        fpLoadFile = open(loadFilePath, 'w')
    except:
        exit('Could not open file for writing %s\n' % loadFilePath)

    try:
        fpQcRpt = open(qcRptPath, 'w')
    except:
        exit('Could not open file for writing %s\n' % qcRptPath)

    # Create lookup of homology IDs to their marker keys
    results = db.sql('''select a.accid, a._object_key as markerKey, m._organism_key
        from acc_accession a, mrk_marker m 
        where a._mgitype_key = 2 
        and a._logicalDB_key in (47, 64, 172)
        and a._object_key = m._marker_key
        and m._marker_status_key = 1''', 'auto')
    print('loading homologyLookup')
    for r in results:
        print('hMrkID: %s orgKey: %s hMrkKey: %s' % (r['accid'], int(r['_organism_key']), int(r['markerKey']) ))
        homologyLookup[r['accid']] = [ int(r['_organism_key']), int(r['markerKey'])]

    # Create lookup of mouse MGI IDs to their marker keys
    results = db.sql('''select a.accid, a._object_key as markerKey
        from acc_accession a, mrk_marker m
        where a._mgitype_key = 2
        and a._logicalDB_key = 1
        and a.prefixPart = 'MGI:'
        and a._object_key = m._marker_key
        and m._marker_status_key = 1''', 'auto')
    for r in results:
        mouseLookup[r['accid']] = r['markerKey']

    return

def process():
    # Purpose: parse file and create load ready file
    # Returns: 0
    # Assumes: lookups have been initialized
    # Effects: Reads file in file system
    # Throws: Nothing

    global rptOne, rptTwo

    # {mgiID:[list of homology IDs], ...}
    homologyDict = {}
    lineCt = 0
    for line in fpInFile.readlines():
        line = str.strip(line)
        lineCt += 1
        #print('lineCt: %s' % lineCt)
        #print('line: %s' % line)
        # skip the header
        if str.find(line, '#') == 0:
            #print(line)
            continue

        tokens = str.split(line)
        mgiID = tokens[0]
        if str.find(mgiID, 'MGI:') != 0:
            continue
        homologyID = tokens[5]
        if not(str.find(homologyID, 'ZFIN:') == 0 or str.find(homologyID, 'HGNC:') == 0 or str.find(homologyID, 'RGD:') == 0):
            continue
        print('\n\nlineCt: %s' % lineCt)
        print('line: %s' % line)
        print('mgiID: %s' % mgiID)
        print('homologyID: %s' % homologyID)
        #Alliance adds a prefix to the zfin id, remove it
        if str.find(homologyID, 'ZFIN:') == 0:
            homologyID = homologyID[5:]

        # Check IDs in the database
        notIn = 0
        mouseKey = ''
        homologyKey = ''
        
        # report if mgiID not in db
        if mgiID not in mouseLookup:
            print('mgiID not in db: %s' % mgiID)
            rptOne = '%s%s%s%s%s' % (rptOne, lineCt, TAB, line, CRT)
            notIn = 1
        else:   
            mouseKey = mouseLookup[mgiID]

        # report if homology ID not in MGI
        if homologyID not in homologyLookup:
            print('homologyID not in db: %s' % homologyID)
            rptTwo = '%s%s%s%s%s' % (rptTwo, lineCt, TAB, line, CRT)
            notIn = 1
        # if either mgiID or homologyID not in MGI skip
        if notIn == 1:
            continue

        # create list of homologies for mouseKey
        if mouseKey not in homologyDict:
            homologyDict[mouseKey] = []
            # add the mouse key, it is second
            homologyDict[mouseKey].append([1, int(mouseKey)])

        l = copy.deepcopy(list(homologyLookup[homologyID])) # l is 2-element list [orgKey, homologyKey]
        print ('mouseKey: %s' % mouseKey)
        print('homologyLookup[homologyID] before: %s' % l)
        idx = organismOrder.index(l[0]) # get the index of the orgKey in the organismOrder list
        print('idx: %s' % idx)
        l[0] = idx                      # and replace the orgKey with the index
        print('l after: %s' % l)
        homologyDict[mouseKey].append(l)# add homologies to the dictionary
        print('homologyLookup[homologyID] after: %s' % homologyLookup[homologyID])

    # now iterate through the clusters and write to the load ready file
    for mKey in homologyDict:
        homologyList = homologyDict[mKey]  # list of lists e.g. [ [2, rKey], [0, hKey], [1, mKey] ]
        print('homologyList: %s' % homologyList)
        sortedList = sorted(homologyList, key=lambda hom: hom[0])    # sort by the index in position 1
        print('sortedList: %s' % sortedList)

        # get the list of marker keys for writing out to load ready file
        keyList = []
        for l in sortedList:
            keyList.append(str(l[1]))
        keyString = ', '.join(keyList)
        print('writing keyString to file: %s' % keyString)
        fpLoadFile.write ('%s%s%s%s' % ('', TAB, keyString, CRT))

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

print('%s' % mgi_utils.date())

print('initializing')
init()

print('processing direct homologies')
process()

print('writing reports')
writeReports()

print('closing files')
closeFiles()

print('%s' % mgi_utils.date())

