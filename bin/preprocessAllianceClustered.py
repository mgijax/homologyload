
##########################################################################
#
# Purpose:
#       From Alliance input file create load ready file
#
# Usage: preprocessAllianceClustered.py
#
# Inputs:
#	1. Alliance file tab-delimited in following format:
#           1. Gene1ID
#           2. Gene1Symbol (not used)
#           3. Gene1SpeciesTaxonID (not used)
#           4. Gene1SpeciesName (not used)
#           5. Gene2ID
#           6. Gene2Symbol (not used)
#           7. Gene2SpeciesTaxonID (not used)
#           8. Gene2SpeciesName (not used)
#           9-13 (not used)
#
#	2. Configuration - see  alliance_clusteredload.config
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
# sc   01/21/2021
#       - initial implementation
###########################################################################

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

# input file from Alliance
inFilePath = os.environ['INPUT_FILE']
print('inFilePath: %s' % inFilePath)

# path to the file to be clustered by the clusterizer
# this was for debugging this first clustered load. I chose to leave this code
# in for future debugging purposes
clustererFilePath = os.environ['INPUT_FILE_CLUSTERER']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

#
# The QC report
qcRptPath = os.environ['QC_RPT']

# Mouse MGI ID/HGNC ID associations from the file
# {mouse MGI ID:[list of human hgncIDs], ...}
mouseMgiToHGNCDict = {}

# HGNC ID/Mouse Marker associations from the database
# {hgncID:marker key, ...}
hgncToMarkerDict = {}

# MGI ID/Mouse Marker associations from the database
# {mgiID:marker key, ...}
mgiToMarkerDict = {}

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'

rptOne = 'Lines where a Mouse MGI ID not in database %s%s%s' % (CRT, sep, CRT)
rptOne = rptOne + 'LineNum%sline%s' % (TAB, CRT)
rptTwo = '%s%sLines where a HGNC ID not in database%s%s%s%s' % (CRT, CRT, CRT, CRT, sep, CRT)
rptTwo = rptTwo + 'LineNum%sline%s' % (TAB, CRT)
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


    global hgncToMarkerDict, mgiToMarkerDict
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


    # get all human markers that have hgncIDs
    results = db.sql('''select distinct a.accid as hgncID, m._Marker_key
        from ACC_Accession a, MRK_Marker m
        where a._MGIType_key = 2
        and a._LogicalDB_key = 64
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        and m._Marker_Status_key = 1
        and m._Organism_key = 2''', 'auto')
    #
    # create hgncID to marker lookup from database
    #
    for r in results:
        hgncID = r['hgncID']
        markerKey = r['_Marker_key']
        hgncToMarkerDict[hgncID] = markerKey

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
    # Purpose: parse file into dictionary 
    # Returns: 0
    # Assumes: mouseMgiToHGNCDict has been initialized
    # Effects: Reads file in file system
    # Throws: Nothing

    global mouseMgiToHGNCDict

    # ignore header lines which start with '#'
    headers = []
    for line in fpInFile.readlines():
        line = str.strip(line)
        if str.find(line, '#') == 0:
            continue
        elif str.find(line, 'Gene1ID') == 0:
            headers = str.split(line, TAB)
            print("headers.index('Gene1ID'): %s" % headers.index('Gene1ID'))
            continue
        tokens =  str.split(line, TAB)
        mgiID = tokens[headers.index('Gene1ID')]
        if str.find(mgiID, 'MGI:') != 0:
            continue
        #print('line: %s' % line)
        #print('mgiID: %s' % mgiID)

        homologyID = tokens[headers.index('Gene2ID')]
        if not str.find(homologyID, 'HGNC:') == 0:
            continue
        #print('homologyID: %s' % homologyID)
        if not mgiID in mouseMgiToHGNCDict:
            mouseMgiToHGNCDict[mgiID] = []
        # add the human  homology to the dictionary
        mouseMgiToHGNCDict[mgiID].append(homologyID)
    return

def process():
    # Purpose: Create load ready file from the Alliance file and the database
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

    for mgiID in mouseMgiToHGNCDict:
        hgncIDList = mouseMgiToHGNCDict[mgiID]
        lineCt += 1
        
        # 1 means error on this line
        error = 0

        clusterFileLine = ''

        # current cluster - if there are no errors it will be added to
        # 'toClusterList'
        currentClusterList = []
        # report and skip lines where hgncID not in the database
        if  mgiID not in list(mgiToMarkerDict.keys()):
            toReport = '%s%s%s%s' % (mgiID, TAB, ''.join(hgncIDList), CRT)
            rptOne = '%s%s%s%s' % (rptOne, lineCt, TAB, toReport)
        
            # if mgiID not in database continue to next input line
            continue

        #
        # mgiID is in the database; check the hgnc IDs
        #

        # add clusters with human to the list
        for id in hgncIDList:
            id = str.strip(id)
            # report and skip lines with hgncId not in the database
            if id not in list(hgncToMarkerDict.keys()): 
                error = 1
                toReport = '%s%s%s%s' % (mgiID, TAB, id, CRT)
                rptTwo = '%s%s%s%s' % (rptTwo, lineCt, TAB, toReport)
                # No need to check any more ids, get out of the loop
                break
            else:
                currentClusterList.append([mgiID, id])
            clusterFileLine = ('%s%s%s%s%s' % \
                (clusterFileLine, mgiID, TAB, id, CRT))

        # if any hgnc IDs not in database continue to next input line
        if error == 1:
                continue

        # if we have a cluster add it to the cluster list and to the file
        if currentClusterList != []:
            # no errors so append the next cluster
            toClusterList = toClusterList + currentClusterList

            fpClustererFile.write(clusterFileLine)
            # if we get here, we know mgiID is in the database and ALL the
            # hgncIds are in the database
                
    fpClustererFile.close()

    clusterDict = clusterize.cluster(toClusterList, 'Alliance')

    # now resolve the ids to database keys; human and mouse gene keys
    for clusterId in list(clusterDict.keys()):
        idTuple = clusterDict[clusterId]
        humanKeyList = []
        mouseKeyList = []
        for id in idTuple:
            if id.startswith('MGI:'):
                mouseKeyList.append(str(mgiToMarkerDict[id]))
            else:
                humanKeyList.append(str(hgncToMarkerDict[id]))
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

print('%s' % mgi_utils.date())

print('initializing')
init()

print('processing clusters')
process()

print('writing reports')
writeReports()

print('closing files')
closeFiles()

print('%s' % mgi_utils.date())
