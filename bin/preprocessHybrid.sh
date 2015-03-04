#!/bin/sh 

#
# Convenience script to use during development to run the preprocessor
#
# This script is a wrapper around the process that determines the hybrid
# clusters and generates the load-ready hybrid input file
#
# Usage:
#
#     preprocessHybrid.sh
#

cd `dirname $0`/..
CONFIG_LOAD=`pwd`/hybridload.config
CONFIG_COMMON=`pwd`/common.config

cd `dirname $0`
LOG=`pwd`/hybridload.log
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
# check that INPUT_FILE_LOAD has been set
#
if [ "${INPUT_FILE_LOAD}" = "" ]
then
    # set STAT for endJobStream.py
    STAT=1
    checkStatus ${STAT} "INPUT_FILE_LOAD not defined"
fi

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
# Run Preprocessor
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo 'Running Preprocessor' >> ${LOG_DIAG}
${HOMOLOGYLOAD}/bin/preprocessHybrid.py
STAT=$?
checkStatus ${STAT} "${HOMOLOGYLOAD}/bin/preprocessHybrid.py"

#
# run postload cleanup and email logs
#
shutDown
