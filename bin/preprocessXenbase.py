
##########################################################################
#
# Purpose:
#       From Xenbase input files create load ready file
#
# Usage: preprocessXenbase.py
#
# Inputs:
#	1. Xenbase expression file tab delimited:
#	    1. Xenbase Gene ID 	
# 	    2+ not used
#
#       2. Xenbase  egID to Xenbase GEne ID file tab delimited:
#           1. Xenbase Gene ID
#           2. Xenopus EG ID
#           3+ not used
#
#	3. Xenbase Gene ID to GenePage ID translation file tab delimited:
# 	    1. Xenbase GenePage ID 	
# 	    2. Xenbase Gene ID
#
#       4. Xenbase Orthology file tab delimited:
#           1. Mouse EG ID
#           2. Xenbase GenePage ID
#
#	5. Configuration - see zfinload.config
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
import Set
import mgi_utils
import clusterize
import db

###--- globals ---###

# constants
TAB= '\t'
CRT = '\n'

# set to 1 for debug to stdout
debug = 0

# EG ID/Xenopus Marker associations from the database
# {egID:marker key, ...}
egToXenMarkerDict = {}

# EG ID/Mouse Marker associations from the database
# {egID:mouse marker key, ...}
egToMouseMarkerDict = {}

#
# paths to input and output files
#

# input files from Xenbase
inFileEgPath = os.environ['INPUT_FILE_EG']
inFileTransPath = os.environ['INPUT_FILE_TRANS']
inFileOrthoPath = os.environ['INPUT_FILE_ORTHO']
# The expression input file
inFileExprPath = os.environ['INPUT_FILE']

# This is the cleaned up load-ready input file
loadFilePath = os.environ['INPUT_FILE_LOAD']

# Xenbase gene IDs from the expression file
exprSet = set([])

# Xenbase gene IDs mapped to their EG IDs from EG file
egDict = {}

# Xenbase gene IDs mapped to their GenePage IDs from the translation file
# {gId: [list of gpIds], ...}
transDict = {}

# Xenbase GenePage IDs mapped to one or more mouse MGI IDs from orthology file
mouseDict = {}

#
# The QC report 
qcRptPath = os.environ['QC_RPT']

# QC report descriptions and column headings
sep = '--------------------------------------------------\n'

rptOne = 'Xenbase Gene IDs from expression file not in the translation file%s%s%s' % (CRT, CRT, sep)

rptTwo = '%s%sXenbase Gene IDs from expression file not in Xenopus EG file%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptThree = '%s%sXenbase Gene Page IDS from the translation file not in the Orthology File%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptFour = '%s%sXenopus EG IDs that map to > 1 Xenbase Gene ID%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptFive = '%s%sXenopus EG IDs not in database%s%s%s' % (CRT, CRT, CRT, CRT, sep)

rptSix = '%s%sMouse EG IDs not in database%s%s%s' % (CRT, CRT,CRT, CRT, sep)

#
# file descriptors
#

fpEgFile = ''
fpTransFile = ''
fpOrthoFile = ''
fpExprFile = ''

# idList replaces this file, but keeping for time being
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

    global egToXenMarkerDict, egToMouseMarkerDict
    global fpEgFile, fpTransFile, fpOrthoFile, fpExprFile
    global fpLoadFile, fpQcRpt, mouseEgMultiGeneIdSet

    mouseEgMultiGeneIdSet = set([])

    user = os.environ['MGD_DBUSER']
    passwordFileName = os.environ['MGD_DBPASSWORDFILE']
    db.useOneConnection(1)
    db.set_sqlUser(user)
    db.set_sqlPasswordFromFile(passwordFileName)

    try:
        fpEgFile = open(inFileEgPath, 'r')
    except:
        exit('Could not open file for reading %s\n' % inFileGenePath)

    try:
        fpTransFile = open(inFileTransPath, 'r')
    except:
        exit('Could not open file for reading %s\n' % inFileTransPath)

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


    # get all xenopus tropicalis markers that are associated with egIds
    results = db.sql('''select distinct a.accid as egId, m._Marker_key
        from ACC_Accession a, MRK_Marker m
        where a._MGIType_key = 2
        and a._LogicalDB_key = 55
        and a.preferred = 1
        and a._Object_key = m._Marker_key
        and m._Marker_Status_key = 1
        and m._Organism_key = 95''', 'auto')
    #
    # create Xenopus egID to marker lookup from database
    #
    for r in results:
        egId = r['egId']
        markerKey = r['_Marker_key']
        egToXenMarkerDict[egId] = markerKey

    # mouse egID to marker lookup from database
    results = db.sql('''select distinct a.accID as egId, m._Marker_key
            from ACC_Accession a, MRK_Marker m
            where a._MGIType_key = 2
            and a._LogicalDB_key = 55
            and a._Object_key = m._Marker_key
            and m._Marker_Status_key = 1
            and m._Organism_key = 1''', 'auto')

    # removed per Richard
    # and a.preferred = 1
    for r in results:
            egId = r['egId']
            markerKey = r['_Marker_key']
            egToMouseMarkerDict[egId] = markerKey

    return

def processInputFiles():
    # Purpose: create data structures from the input files
    # Returns: 0
    # Assumes: Nothing
    # Effects: None
    # Throws: Nothing

    global exprSet, egDict, mouseDict, transDict
    global xenEgToGeneIdDict

    xenEgToGeneIdDict = {}

    for line in fpExprFile.readlines():
        tokens = string.split(line, TAB)
        # get the xbgid
        exprSet.add(string.strip(tokens[0]))

    #
    # Xenopus tropicalis gene ID to EG ID file
    #
    for line in fpEgFile.readlines():
        tokens = string.split(line, TAB)
        # get the Xenbase Gene ID
        gId = string.strip(tokens[0])
        # get the Xenopus EG ID
        egId = string.strip(tokens[2])
        if egId == '':
            egId = 'None'
        egDict[gId] = egId
        if egId !='None':
            if not xenEgToGeneIdDict.has_key(egId): 
                xenEgToGeneIdDict[egId] = []
            xenEgToGeneIdDict[egId].append(gId)
    for egId in xenEgToGeneIdDict.keys():
        if len(xenEgToGeneIdDict[egId]) > 1:
            # write to bad egId report
            geneIds = xenEgToGeneIdDict[egId]
            mouseEgMultiGeneIdSet.add(egId)
            # now remove the gene ID from egDict because it participates in an
            # eg ID that maps to multiple gene IDs
            for g in geneIds:
                egDict.pop(g)
    if debug:
        print 'Xenopus tropicalis gene ID to EG ID file'
        keys = sorted(egDict.keys())
        for key in keys:
            print 'gId: %s egId: %s' % (key, egDict[key])
    #
    # Xenbase Gene Page ID to list of Xenbase Gene IDs file
    #
    for line in fpTransFile.readlines():
        tokens = string.split(line, TAB)
        # get the Xenbase Gene Page ID
        gpId = string.strip(tokens[0])
        # get the Xenbase Gene IDs from every other field, ignoring
        # the gene symbols
        for gId in tokens[2::2]:
            transDict[gId] = gpId
    if debug:
        print 'Xenbase Gene ID to Gene Page ID translation file'
        keys = sorted(transDict.keys())
        for key in keys:
            print 'gId: %s gpId: %s' % (key, transDict[key])
    #
    # Xenopus Gene Page ID to mouse EG ID  file
    #
    for line in fpOrthoFile.readlines():
        tokens = string.split(line, TAB)
        # get the mouse EG ID
        egID = string.strip(tokens[0])
        # get the Xenbase Gene Page ID
        gpId = string.strip(tokens[1])
        # one mouse egId to many genePage IDs
        mouseDict[gpId] = egID
    if debug:
        print 'Xenopus Gene Page ID to mouse EG ID file'
        keys = sorted(mouseDict.keys())
        for key in keys:
            print 'gpId: %s mouseEgId: %s' % (key, mouseDict[key])

def process():
    # Purpose: Create load ready file and  QC reports from Xenbase files 
    #	and the database
    # Returns: 0
    # Assumes: All lookup structures have been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    global noTransSet, noXenEgSet, noOrthSet, mouseNotInDbSet, xenNotInDbSet
    global rptOne, rptTwo, rptThree, rptFour, rptFive, rptSix
    # database lookups
    # egToXenMarkerDict, egToMouseMarkerDict
    # input file lookups
    # egDict = {}
    # transDict = {}
    # mouseDict = {}

    # dictionary of id pairs to send to the clusterizer
    toClusterList = []

    # set of xen geneIds with no translation
    noTransSet = set([])

    # set of xen geneIds with no xen eg id 
    noXenEgSet = set([])
 
    # set of geneIds that map to trans file, but gpId doesn't map to orth file
    noOrthSet = set([])

    # set of xen egIds not in the database
    xenNotInDbSet = set([])

    # set of mouse egIds not in the database
    mouseNotInDbSet = set([])

    for geneId in exprSet:
        # Join to trans and eg file on geneId
        if not  transDict.has_key(geneId):
            noTransSet.add(geneId)
            continue
        gpId = transDict[geneId]
        skip = 0
        if not egDict.has_key(geneId):
            noXenEgSet.add(geneId)
            skip = 1
        elif egDict[geneId] == 'None':
            noXenEgSet.add('%s: %s' % (geneId, egDict[geneId]))
            skip = 1
        if skip == 1:
            continue
        xenEg = egDict[geneId]
        # join from trans to orth file on genePageId
        if not mouseDict.has_key(gpId):
            noOrthSet.add(gpId)
            continue
        mouseEg = mouseDict[gpId]
        notInDb = 0
        if not egToXenMarkerDict.has_key(xenEg):
            xenNotInDbSet.add(xenEg)
            notInDb = 1
        if not egToMouseMarkerDict.has_key(mouseEg):
            mouseNotInDbSet.add(mouseEg)
            notInDb = 1
        if notInDb:
            continue
        # now we have a homology
        toClusterList.append([xenEg, mouseEg])

    # send cluster list to clusterizer
    clusterDict = clusterize.cluster(toClusterList, 'XENBASE')
    # now resolve the ids to database keys; xenbase and mouse gene keys
    # NOT SORTING BY ORGANISM for sequenceNum because we don't need to.
    # If we find we need for this load we will need to create a mouse and a
    # xenopus lookup by EG ID to determine which organism in order to
    # order correctly (mouse first then xenopus)
    for clusterId in clusterDict.keys():
        idTuple = clusterDict[clusterId]
        mouseKeyList = []
        xenKeyList = []
        for id in idTuple:
            if id in egToMouseMarkerDict.keys():
                mouseKeyList.append(str(egToMouseMarkerDict[id]))
            elif id in egToXenMarkerDict.keys():
                xenKeyList.append(str(egToXenMarkerDict[id]))
            else:
                print 'not xenopus or mouse'
        keyList = mouseKeyList + xenKeyList
        keyString = ', '.join(keyList)
        fpLoadFile.write('%s%s%s%s' % (clusterId, TAB, keyString, CRT))

    for id in noTransSet:
        rptOne = '%s%s%s' % (rptOne, id, CRT)
    rptOne = '%s%sTotal IDs: %s%s' % (rptOne, CRT, len(noTransSet), CRT)

    for id in noXenEgSet:
        rptTwo = '%s%s%s' % (rptTwo, id, CRT)
    rptTwo = '%s%sTotal IDs: %s%s' % (rptTwo, CRT, len(noXenEgSet), CRT)

    for id in noOrthSet:
        rptThree = '%s%s%s' % (rptThree, id, CRT)
    rptThree = '%s%sTotal IDs: %s%s' % (rptThree, CRT, len(noOrthSet), CRT)

    for id in mouseEgMultiGeneIdSet:
        rptFour = '%s%s%s' % (rptFour, id, CRT)
    rptFour = '%s%sTotal IDs: %s%s' % (rptFour, CRT, len(mouseEgMultiGeneIdSet), CRT)

    for id in mouseNotInDbSet:
        rptFive = '%s%s%s' % (rptFive, id, CRT)
    rptFive = '%s%sTotal IDs: %s%s' % (rptFive, CRT, len(mouseNotInDbSet), CRT)

    for id in xenNotInDbSet:
        rptSix = '%s%s%s' % (rptSix, id, CRT)
    rptSix = '%s%sTotal IDs: %s%s' % (rptSix, CRT, len(xenNotInDbSet), CRT)

    return

def writeReports():
    # Purpose: writes out all sections of the QC report
    # Returns: 0
    # Assumes: rptOne - rptSix has been initialized
    # Effects: Writes to the file system
    # Throws: Nothing

    fpQcRpt.write(rptOne)
    fpQcRpt.write(rptTwo)
    fpQcRpt.write(rptThree)
    fpQcRpt.write(rptFour)
    fpQcRpt.write(rptFive)
    fpQcRpt.write(rptSix)

    return

def closeFiles():
    # Purpose: closes file descriptors and database connection
    # Returns: 0
    # Assumes: file descriptors have been initialized
    # Effects:  None
    # Throws: Nothing

    fpEgFile.close()
    fpOrthoFile.close()
    fpExprFile.close()
    fpLoadFile.close()
    fpQcRpt.close()

    # close the database connection
    db.useOneConnection(0)
    
    return

###--- main program ---###

print '%s' % mgi_utils.date()

print 'initializing'
init()

print 'processing input files'
processInputFiles()

print 'processing clusters'
process()

print 'writing reports'
writeReports()

print 'closing files'
closeFiles()

print '%s' % mgi_utils.date()
