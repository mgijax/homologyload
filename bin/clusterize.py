#!/usr/local/bin/python

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
	idOne = pair[0]
	idTwo = pair[1]
	if not idOneToIdTwoDict.has_key(idOne):
	    idOneToIdTwoDict[idOne] = []
	idOneToIdTwoDict[idOne].append(idTwo)
	if idTwo != 'None':
	    if not idTwoToIdOneDict.has_key(idTwo):
		idTwoToIdOneDict[idTwo] = []
	    idTwoToIdOneDict[idTwo].append(idOne)
    # set of all clusters (uniq set of tuples)
    clusterSet = set()

    # current number of clusters created, for generating cluster ID
    clusterCt = 0

    # iterate through column 1 ids
    for idOne in idOneToIdTwoDict.keys():
	# current cluster
	currentList = []
	# this includes input rows where idOne exists, but no idTwo
	currentList.append(idOne)
	# add all column 2 IDs to the cluster
	for idTwo in idOneToIdTwoDict[idOne]:
	    if idTwo != 'None' and idTwo not in currentList:
		currentList.append(idTwo)
	    # add all column 1 IDs to the cluster
	    if idTwoToIdOneDict.has_key(idTwo):
		for hId in idTwoToIdOneDict[idTwo]:
		    if hId not in currentList:
			currentList.append(hId)
		    # check idOneToTwoDict for hId and add it's values if found and not already in the list
		    if hId in idOneToIdTwoDict:
			iList = idOneToIdTwoDict[hId]
			for i in iList:
			    if i not in currentList:
				currentList.append(i)
	# sort the list so dups will not be created in clusterSet
	currentList.sort()

	# add list converted to tuple to clusterSet
	clusterSet.add(tuple(currentList))
    # now assign ids to the clusters creating a dict to return
    # we'll remove IDs when we change the schema and simply return a list
    # of tuples
    # {clusterID:(id1, ..., idn), ...} 
    clusterDict = {}
    clusterCt = 0
    for c in clusterSet:
	clusterCt += 1
	nextId = '%s:%s' % (idPrefix, clusterCt)
	#clusterDict[nextId] = ', '.join(c)
	clusterDict[nextId] = c
	
    return clusterDict
