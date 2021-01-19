#
#  checkColumns.py
###########################################################################
#
#  Purpose:
#
#       This script checks that there are the correct number of columns
#	in a file
#
#  Usage:
#
#      checkColumns.py  filename numColumns	
#
#      where:
#          filename = path to the input file
#
#  Env Vars:
#
#      The following environment variables are set by the configuration
#      files that are sourced by the wrapper script:
#
#  Inputs:
#
#  Outputs:
#
#  Exit Codes:
#
#      0:  Successful completion
#      1:  An exception occurred
#
#  Implementation:
#
#  Notes:  None
#
###########################################################################

import string
import sys

USAGE = 'Usage: checkColumns.py  inputFile numColumns'
TAB = '\t'

inputFile = None
fpInput = None
numColumns = None
errors = 0

def checkArgs ():
    # Purpose: Validate the arguments to the script.
    # Returns: 1 if input file has <> 3 columns
    # Assumes: Nothing
    # Effects: Sets global variables.
    # Throws: Nothing

    global inputFile, numColumns
    if len(sys.argv) != 3:
        print(USAGE)
        sys.exit(1)

    inputFile = sys.argv[1]
    numColumns = int(sys.argv[2])
    return

def openFile ():
    # Purpose: Open the file for reading
    # Returns: 1 if file descriptor cannot be initialized
    # Assumes: Nothing
    # Effects: Sets global variables.
    # Throws: Nothing

    global fpInput

    try:
        fpInput = open(inputFile, 'r')
    except:
        print('Cannot open input file: ' + inputFile)
        sys.exit(1)
    return

def checkColumns ():
    # Purpose: check the file for proper number of columns
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Sets global variables.
    # Throws: Nothing

    global errors
    lineNum = 0
    for line in fpInput.readlines():
        if str.find(line, '#') == 0:
            continue
        colError = 0
        lineNum = lineNum + 1
        columns = str.split(line, TAB)
        # remove newline from last column
        last = columns[-1].strip()
        columns[-1] = last
        nc = len(columns) 
        if nc < numColumns:
            print('Missing Column(s) in %s on line %s: %s ' % (inputFile, lineNum, columns))
            continue
    return

def closeFile():
    # Purpose: Close the files.
    # Returns: Nothing
    # Assumes: Nothing
    # Effects: Nothing
    # Throws: Nothing

    global fpInput

    fpInput.close()
    return

checkArgs()
openFile()
checkColumns()
closeFile()
if errors > 0:
    sys.exit(1)
sys.exit(0)
