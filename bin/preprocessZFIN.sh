#!/bin/sh

#
# Convenience script to use during development to run the preprocessor
#
# This script is a wrapper around the process that QC's the ZFIN
# input file and generates a load-ready input file
#
# Usage:
#
#     preprocessZFIN.sh
#

#
# Create a temporary file and make sure that it is removed when this script
# terminates.
#
TMP2_FILE=/tmp/`basename $0`.$$
touch ${TMP2_FILE}
trap "rm -f ${TMP2_FILE}" 0 1 2 15


#
# check that two additional default input files  have been set
# Note that the expression file is the default input file and is
# checked in the calling wrapper script
# 
if [ "${INPUT_FILE_GENE_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_GENE_DEFAULT not defined"
fi

if [ "${INPUT_FILE_ORTHO_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_ORTHO_DEFAULT not defined"
fi

#
# check that two additional input files have been set
# 
if [ "${INPUT_FILE_GENE}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_GENE not defined"
fi

if [ "${INPUT_FILE_ORTHO}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_ORTHO not defined"
fi

# copy the additional two files from /data/downloads to the input dir
echo "copying ${INPUT_FILE_GENE_DEFAULT} to ${INPUTDIR}"
cp -p ${INPUT_FILE_GENE_DEFAULT} ${INPUTDIR}
echo "copying ${INPUT_FILE_ORTHO_DEFAULT} to ${INPUTDIR}"
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
echo "Running Sanity Checks on additional two input files" >> ${LOG_DIAG}
FILE_ERROR=0

echo '                         Sanity Errors in additional two input files' >> ${SANITY_RPT}
echo '---------------------------------------------------------------' >> ${SANITY_RPT}
echo ''

#
# check each of additional two input files for length, remove whitespace
#
len=`cat ${INPUT_FILE_GENE} | wc -l | sed 's/ //g'`
if [ ${len} -lt ${MIN_LENGTH} ]
then
   echo "\n\nInput file ${INPUT_FILE_GENE} does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
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

# gene.txt - column 1: zfin gene id, column 4: NCBI gene id
COLUMNS=4
echo "checkColumns ${INPUT_FILE_GENE} ${COLUMNS}"
checkColumns ${INPUT_FILE_GENE} ${COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

# mouse_orthos.txt - column 1: zfin gene id, column 5: MGI mouse gene id
COLUMNS=5
echo "checkColumns ${INPUT_FILE_ORTHO} ${COLUMNS}"
checkColumns ${INPUT_FILE_ORTHO} ${COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

#
# If the input file had sanity errors exit
#
echo "TMP2_FILE: ${TMP2_FILE}"
echo "Contents: "
cat ${TMP2_FILE}
if [ ${FILE_ERROR} -ne 0 ]
then
    cat ${TMP2_FILE} >> ${SANITY_RPT}
    echo "Sanity errors in input file. See ${SANITY_RPT}" >> ${LOG_DIAG}
    echo "Sanity errors in input file. See ${SANITY_RPT}" >> ${LOG_PROC}
    
    exit 1
else
    echo "No sanity errors in input file" >> ${SANITY_RPT}
fi

#
# QC checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running QC Checks' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/preprocessZFIN.py
STAT=$?
exit ${STAT}

exit 0
