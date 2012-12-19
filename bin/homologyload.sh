#!/bin/sh -x

#
# load Homologene
#
# Usage:
#
#     homologyload.sh
#

cd `dirname $0`/..
CONFIG_LOAD=`pwd`/homologyload.config

cd `dirname $0`
LOG=`pwd`/homologyload.log
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
# verify & source the configuration file
#

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
    ${HOMOLOGYLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS}  ${TMP_FILE}
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
# There should be a "lastrun" file in the input directory that was created
# the last time the load was run for this input file. If this file exists
# and is more recent than the input file, the load does not need to be run.
#
LASTRUN_FILE=${INPUTDIR}/lastrun
if [ -f ${LASTRUN_FILE} ]
then
    if /usr/local/bin/test ${LASTRUN_FILE} -nt ${INPUT_FILE_DEFAULT}
    then

        echo "Input file has not been updated - skipping load" | tee -a ${LOG_PROC}
        # set STAT for shutdown
        STAT=0
        echo 'shutting down'
        shutDown
        exit 0
    fi
fi

# 
# Sanity checks
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
#echo "'$len'"
if [ ${len} -lt ${MIN_LENGTH} ] 
then
   echo "\n\nInput file does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
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

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running QC Checks' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/runQC.py
STAT=$?
checkStatus ${STAT} "${HOMOLOGYLOAD}/bin/runQC.py"

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running homologyload.py' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/homologyload.py
STAT=$?
checkStatus ${STAT} "${HOMOLOGYLOAD}/bin/homologyload.py"

COLDELIM="\t"
LINEDELIM="\n"

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'BCP data into MRK_Cluster'  >> ${LOG_DIAG}

TABLE=MRK_Cluster

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object >> ${LOG_DIAG}

# BCP new data 
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object >> ${LOG_DIAG}

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'BCP data into MRK_ClusterMember'  >> ${LOG_DIAG}

TABLE=MRK_ClusterMember

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object >> ${LOG_DIAG}

# BCP new data 
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object >> ${LOG_DIAG}

echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'BCP data into ACC_Accession'  >> ${LOG_DIAG}

TABLE=ACC_Accession

# Drop indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_drop.object >> ${LOG_DIAG}

# BCP new data
${MGI_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}

# Create indexes
${MGD_DBSCHEMADIR}/index/${TABLE}_create.object >> ${LOG_DIAG}

#
# Touch the "lastrun" file to note when the load was run.
#
#if [ ${STAT} = 0 ]
#then
#    touch ${LASTRUN_FILE}
#fi

#
# run postload cleanup and email logs
#
shutDown

