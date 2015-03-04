#!/bin/sh

#
# Wrapper to run the Xenbase preprocessor. It copies the three additional input 
# files and does sanity checks on them, then runs the python preprocessor
#
# Usage:
#
#     preprocessXenbase.sh
#

#
# Create a temporary file and make sure that it is removed when this script
# terminates.
#
TMP2_FILE=/tmp/`basename $0`.$$
touch ${TMP2_FILE}
trap "rm -f ${TMP2_FILE}" 0 1 2 15

#
# check that three additional default input files  have been set
# Note that the expression file is the default input file and is
# checked in the calling wrapper script
# 
if [ "${INPUT_FILE_EG_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_EG_DEFAULT not defined"
fi

if [ "${INPUT_FILE_TRANS_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_TRANS_DEFAULT not defined"
fi

if [ "${INPUT_FILE_ORTHO_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_ORTHO_DEFAULT not defined"
fi

#
# check that three additional input files have been set
# 
if [ "${INPUT_FILE_EG}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_EG not defined"
fi

if [ "${INPUT_FILE_TRANS}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_TRANS not defined"
fi

if [ "${INPUT_FILE_ORTHO}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_ORTHO not defined"
fi

# copy the additional three files from /data/downloads to the input dir
cp -p ${INPUT_FILE_EG_DEFAULT} ${INPUTDIR}
cp -p ${INPUT_FILE_TRANS_DEFAULT} ${INPUTDIR}
cp -p ${INPUT_FILE_ORTHO_DEFAULT} ${INPUTDIR}

#
# FUNCTION: Check for lines with missing columns and data in input file and
#           write the line numbers to the sanity report.
#
checkColumns ()
{
    FILE=$1         # The input file to check
    NUM_COLUMNS=$2  # The number of columns expected in each input record
    ${HOMOLOGYLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS} >>  ${TMP2_FILE}
    if [ `cat ${TMP2_FILE} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi
}

#####################################
#
# Main
#
#####################################

#
# Run sanity checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Running Sanity Checks on additional three input files" >> ${LOG_DIAG}
FILE_ERROR=0

echo '                         Sanity Errors in additional three input files' >> ${SANITY_RPT}
echo '---------------------------------------------------------------' >> ${SANITY_RPT}
echo ''

#
# check each of additional three input files for length, remove whitespace
#
len=`cat ${INPUT_FILE_EG} | wc -l | sed 's/ //g'`
if [ ${len} -lt ${MIN_LENGTH} ]
then
   echo "\n\nInput file ${INPUT_FILE_EG} does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
   FILE_ERROR=1
fi

#
# check each of additional three input files for length, remove whitespace
#
len=`cat ${INPUT_FILE_TRANS} | wc -l | sed 's/ //g'`
if [ ${len} -lt ${MIN_LENGTH} ]
then
   echo "\n\nInput file ${INPUT_FILE_TRANS} does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
   FILE_ERROR=1
fi

len=`cat ${INPUT_FILE_ORTHO} | wc -l | sed 's/ //g'`
if [ ${len} -lt ${MIN_LENGTH} ]
then
   echo "\n\nInput file ${INPUT_FILE_ORTHO} does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
   FILE_ERROR=1
fi

#
# check each input file for proper number of columns
#

# column 1: xenbase gene id, column 3: xenopus EG id
COLUMNS=3
checkColumns ${INPUT_FILE_EG} ${COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

# column 1: xenbase gene page id, column 2: xenbase gene id
COLUMNS=2
checkColumns ${INPUT_FILE_TRANS} ${COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

#  column 1: mouse EG id, column 2: Xenbase gene page id
COLUMNS=2
checkColumns ${INPUT_FILE_ORTHO} ${COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

#
# If the input file had sanity errors exit
#
cat ${TMP2_FILE}
if [ ${FILE_ERROR} -ne 0 ]
then
    cat ${TMP2_FILE} >> ${SANITY_RPT}
    echo "Sanity errors in input file. See ${SANITY_RPT}" >> ${LOG_DIAG}
    echo "Sanity errors in input file. See ${SANITY_RPT}" >> ${LOG_PROC}
    
    exit 1
else
    echo "No sanity errors in additional input files" >> ${SANITY_RPT}
fi

#
# QC checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running Preprocessor' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/preprocessXenbase.py
STAT=$?
exit ${STAT}

exit 0
