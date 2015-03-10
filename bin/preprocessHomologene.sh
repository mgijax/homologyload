#!/bin/sh 

#
# Convenience script to use during development to run the preprocessor
#
# This script is a wrapper around the process that QC's the HomoloGene
# input file and generates a load-ready input file
# Usage:
#
#     preprocessHomologene.sh
#

cd `dirname $0`/..
CONFIG_LOAD=`pwd`/homologeneload.config
CONFIG_COMMON=`pwd`/common.config

cd `dirname $0`
LOG=`pwd`/homologeneload.log
rm -rf ${LOG}

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 0 ]
then
    echo ${Usage} | tee -a ${LOG}
    exit 1
fi

#
# Create a temporary file and make sure that it is removed when this script
# terminates.
#
TMP_FILE=/tmp/`basename $0`.$$
touch ${TMP_FILE}
trap "rm -f ${TMP_FILE}" 0 1 2 15

#
# verify & source the configuration files
#

if [ ! -r ${CONFIG_COMMON} ]
then
    echo "Missing configuration file: ${CONFIG_COMMON}"
    exit 1
fi

. ${CONFIG_COMMON}

if [ ! -r ${CONFIG_LOAD} ]
then
    echo "Cannot read configuration file: ${CONFIG_LOAD}"
    exit 1
fi

. ${CONFIG_LOAD}

#
#  Source the DLA library functions.
#

if [ "${DLAJOBSTREAMFUNC}" != "" ]
then
    if [ -r ${DLAJOBSTREAMFUNC} ]
    then
        . ${DLAJOBSTREAMFUNC}
    else
        echo "Cannot source DLA functions script: ${DLAJOBSTREAMFUNC}" | tee -a ${LOG}
        exit 1
    fi
else
    echo "Environment variable DLAJOBSTREAMFUNC has not been defined." | tee -a ${LOG}
    exit 1
fi

#
# check that INPUT_FILE_DEFAULT has been set
#
if [ "${INPUT_FILE_DEFAULT}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_DEFAULT not defined"
fi

#
# check that INPUT_FILE has been set
#
if [ "${INPUT_FILE}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE not defined"
fi

#
# check that INPUT_FILE_LOAD has been set
#
if [ "${INPUT_FILE_LOAD}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_LOAD not defined"
fi

# copy the latest file from /data/downloads to the input dir
cp -p ${INPUT_FILE_DEFAULT} ${INPUTDIR}

#
# FUNCTION: Check for lines with missing columns and data in input file and
#           write the line numbers to the sanity report.
#
checkColumns ()
{
    FILE=$1         # The input file to check
    REPORT=$2       # The sanity report to write to
    NUM_COLUMNS=$3  # The number of columns expected in each input record
    ${HOMOLOGYLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS} >> ${TMP_FILE}
    cat ${TMP_FILE} | tee -a ${REPORT}
    if [ `cat ${TMP_FILE} | wc -l` -eq 0 ]
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
# createArchive including OUTPUTDIR, startLog, getConfigEnv
# sets "JOBKEY"
preload ${OUTPUTDIR}

#
# rm all files/dirs from OUTPUTDIR
#
cleanDir ${OUTPUTDIR}

#
# Run sanity checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Running Sanity Checks" >> ${LOG_DIAG}
FILE_ERROR=0

echo '                         Sanity Errors' > ${SANITY_RPT}
echo '---------------------------------------------------------------' >> ${SANITY_RPT}
echo ''

checkColumns ${INPUT_FILE} ${SANITY_RPT} ${NUM_COLUMNS}
if [ $? -ne 0 ]
then
    FILE_ERROR=1
fi

# check file length, remove whitespace
len=`cat ${INPUT_FILE} | wc -l | sed 's/ //g'`

if [ ${len} -lt ${MIN_LENGTH} ] 
then
   echo "" >> ${SANITY_RPT}
   echo "" >> ${SANITY_RPT}
   echo "Input file does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
   FILE_ERROR=1
fi

#
# If the input file had sanity errors exit
#
STAT=0
if [ ${FILE_ERROR} -ne 0 ]
then
    echo "Sanity errors in input file. See ${SANITY_RPT}"
    echo "Sanity errors in input file. See ${SANITY_RPT}" >>  ${LOG_DIAG}  ${LOG_PROC}
    # set STAT for shutdown
    STAT=${FILE_ERROR}
    shutDown
    exit 1
fi

#
# QC checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running Preprocessor' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/preprocessHGNC.py
STAT=$?
checkStatus ${STAT} "${HOMOLOGYLOAD}/bin/preprocessHGNC.py"

#
# run postload cleanup and email logs
#
shutDown
