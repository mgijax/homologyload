#!/usr/local/bin/python

##########################################################################
#
# Purpose:
#       Create cluster file with auto-generated IDs
#
# Usage: clusterize.py
#
# Inputs:
#       1. Two column tab-delimited file of ids to be clustered
#       2. Each column containing one ID
#
# Outputs:
#        1. Two column tab-delimited file of clusters
#	 1.1. generated cluster ID
#	 1.2  comma delimited list of cluster member IDs
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

import string
import Set

#def cluster(inputFile, idPrefix):
def cluster(toClusterList, idPrefix):
    #fpCF = open(inputFile, 'r')
    
    # read file into two data structures representing column 1 and column 2
    # e.g. idOne=human NCBI ID, idTwo=mouse MGI ID

    # ids from column 1 mapped to ids in column 2
    idOneToIdTwoDict = {}

    # ids from column 2 mapped to ids from column 2
    idTwoToIdOneDict = {}

    #for line in fpCF.readlines():
    for pair in toClusterList:
	#idOne, idTwo = string.split(line)
	idOne = pair[0]
	idTwo = pair[1]
	if not idOneToIdTwoDict.has_key(idOne):
	    idOneToIdTwoDict[idOne] = []
	idOneToIdTwoDict[idOne].append(idTwo)
	if idTwo != 'None':
	    if not idTwoToIdOneDict.has_key(idTwo):
		idTwoToIdOneDict[idTwo] = []
	    idTwoToIdOneDict[idTwo].append(idOne)
    #fpCF.close()
    print idOneToIdTwoDict
    print idTwoToIdOneDict
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
	    if idTwo != 'None':
		currentList.append(idTwo)
	    # add all column 1 IDs to the cluster
	    if idTwoToIdOneDict.has_key(idTwo):
		for hId in idTwoToIdOneDict[idTwo]:
		    if hId not in currentList:
			currentList.append(hId)

	# sort the list so dups will not be created in clusterSet
	currentList.sort()

	# add list converted to tuple to clusterSet
	clusterSet.add(tuple(currentList))

    # now assign ids to the clusters creating a dict to return
    # we'll remove IDs when we change the schema and simply return a list
    # of lists
    clusterDict = {}
    clusterCt = 0
    for c in clusterSet:
	clusterCt += 1
	nextId = '%s:%s' % (idPrefix, clusterCt)
	clusterDict[nextId] = ', '.join(c)

    return clusterDict
