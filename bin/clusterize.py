
##########################################################################
#
# Purpose:
#       Create cluster file with auto-generated IDs
#
# Usage: clusterize.py
#
# Inputs:
#       1. List of lists
#           [ [idType1, idType2], [idType1, idType2], ...]
#
# Outputs:
#        1. Dictionary of auto-generated clusterIDs mapped to tuples of clusters
#	     {clusterID:(id1, ..., idn), ...}
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

import string
import Set

###--- functions ---###

# getRoots from:
# http://stackoverflow.com/questions/10301000/python-connected-components
def getRoots(aNeigh):
    def findRoot(aNode,aRoot):
        while aNode != aRoot[aNode][0]:
            aNode = aRoot[aNode][0]
        return (aNode,aRoot[aNode][1])
    myRoot = {}
    for myNode in aNeigh.keys():
        myRoot[myNode] = (myNode,0)
    for myI in aNeigh:
        for myJ in aNeigh[myI]:
            (myRoot_myI,myDepthMyI) = findRoot(myI,myRoot)
            (myRoot_myJ,myDepthMyJ) = findRoot(myJ,myRoot)
            if myRoot_myI != myRoot_myJ:
                myMin = myRoot_myI
                myMax = myRoot_myJ
                if  myDepthMyI > myDepthMyJ:
                    myMin = myRoot_myJ
                    myMax = myRoot_myI
                myRoot[myMax] = (myMax,max(myRoot[myMin][1]+1,myRoot[myMax][1]))
                myRoot[myMin] = (myRoot[myMax][0],-1)
    myToRet = {}
    for myI in aNeigh:
        if myRoot[myI][0] == myI:
            myToRet[myI] = []
    for myI in aNeigh:
        myToRet[findRoot(myI,myRoot)[0]].append(myI)
    return myToRet

def cluster(toClusterList, idPrefix):
    # Purpose: finds clusters of related elements in toClusterList, 
    # a list of two-element lists
    # Returns: A dictionary mapping an auto-generated ID (using idPrefix)
    # to a cluster (tuple) of IDs
    # Assumes: toClusterList is a list of two-element lists
    # Effects: Nothing
    # Throws: Nothing

    #fpCF = open(inputFile, 'r')
    
    # read file into two data structures representing column 1 and column 2
    # e.g. idOne=human NCBI ID, idTwo=mouse MGI ID

    # ids from column 1 mapped to ids in column 2
    idOneToIdTwoDict = {}

    # ids from column 2 mapped to ids from column 2
    idTwoToIdOneDict = {}

    #for line in fpCF.readlines():
    for pair in toClusterList:
        idOne = pair[0] # egID
        idTwo = pair[1] # MGI ID
        #
        if not idOneToIdTwoDict.has_key(idOne):
            idOneToIdTwoDict[idOne] = []
        if idTwo != 'None':
            idOneToIdTwoDict[idOne].append(idTwo)
            if not idTwoToIdOneDict.has_key(idTwo):
                idTwoToIdOneDict[idTwo] = []
            idTwoToIdOneDict[idTwo].append(idOne)
                 
    # concatenate the dicts together
    allDict = dict(idOneToIdTwoDict.items() + idTwoToIdOneDict.items())
    #print 'allDict:'
    #for key in allDict:
        #print 'key: %s value: %s' % (key, allDict[key])

    #print getRoots(allDict)
    clusterDict = getRoots(allDict)
    #print 'clusterDict'
    #for key in clusterDict:
        #print 'key: %s value: %s' % (key, clusterDict[key])
    #print clusterDict
    # {clusterID:(id1, ..., idn), ...} 
    namedDict = {}
    clusterCt = 0
    for c in clusterDict:
        cluster = clusterDict[c]
        
        clusterCt += 1
        nextId = '%s:%s' % (idPrefix, clusterCt)
        #clusterDict[nextId] = ', '.join(c)
        namedDict[nextId] = cluster
        
    return namedDict
