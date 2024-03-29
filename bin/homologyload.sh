#!/bin/sh

#
# This script is a wrapper around the process that QC's a homology
# input file, generates a load-ready input file, and runs the load
#
# Usage:
#
#     homologyload.sh configFile
#
USAGE='Usage: homologyload.sh configFile'

cd `dirname $0`/..

CONFIG_COMMON=`pwd`/common.config
LOG=`pwd`/homologyload.log
rm -rf ${LOG}

#
#  Verify the argument(s) to the shell script.
#
if [ $# -ne 1 ]
then
    echo ${USAGE}; exit 1   
else
    CONFIG_LOAD=$1   
fi
CONFIG_LOAD=`pwd`/${CONFIG_LOAD}

#
# Create a temporary file and make sure that it is removed when this script
# terminates.
#
TMP_FILE=/tmp/`basename $0`.$$
touch ${TMP_FILE}
trap "rm -f ${TMP_FILE}" 0 1 2 15

#
# BCP delimiters
#
COLDELIM="\t"
LINEDELIM="\n"

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
# get the homologene version
#
if [ "${RELEASE_NO_FILE}" != "" ]
then
    HOMOLOGY_VERSION=`cat ${RELEASE_NO_FILE}`

    export HOMOLOGY_VERSION
fi


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

#
# copy file from default to input - some files are gzipped
# loads with no input file will be specified as 'None'
#

if [ "${INPUT_FILE_DEFAULT}" != "None" ]
then

    echo "copying ${INPUT_FILE_DEFAULT} to ${INPUTDIR}" >> ${LOG_DIAG}
    # copy the latest file from /data/downloads to the input dir
    if [ ${INPUT_FILE_DEFAULT: -2} = 'gz' ]
    then
        zcat ${INPUT_FILE_DEFAULT} > ${INPUT_FILE}
    else
        cp ${INPUT_FILE_DEFAULT} ${INPUTDIR}
    fi
fi

#
# FUNCTION: Check for lines with missing columns and data in input file and
#           write the line numbers to the sanity report.
#
checkColumns ()
{
    FILE=$1         # The input file to check
    REPORT=$2       # The sanity report to write to
    NUM_COLUMNS=$3  # The number of columns expected in each input record
    ${PYTHON} ${HOMOLOGYLOAD}/bin/checkColumns.py ${FILE} ${NUM_COLUMNS}
    cat ${TMP_FILE} | tee -a ${REPORT}
    if [ `cat ${TMP_FILE} | wc -l` -eq 0 ]
    then
        return 0
    else
        return 1
    fi
}

#
# FUNCTION: run sanity checks on input file
#           write the line numbers to the sanity report.
#
runSanityChecks()
{
    echo "" >> ${LOG_DIAG}
    date >> ${LOG_DIAG}
    echo "Running Sanity Checks" >> ${LOG_DIAG}
    FILE_ERROR=0

    echo '                         Sanity Errors in Primary Input File' > ${SANITY_RPT}
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
       echo "Input file ${INPUT_FILE} does not have minimum length. Required: ${MIN_LENGTH} Found: ${len}" | tee -a ${SANITY_RPT}
       FILE_ERROR=1
    fi

    #
    # If the input file had sanity errors exit
    #
    STAT=0
    if [ ${FILE_ERROR} -ne 0 ]
    then
	echo "Sanity errors in ${INPUT_FILE}. See ${SANITY_RPT}"
	echo "Sanity errors in ${INPUT_FILE}. See ${SANITY_RPT}" >>  ${LOG_DIAG}  ${LOG_PROC}
	# set STAT for shutdown
	STAT=${FILE_ERROR}
	shutDown
	exit 1
    else
	echo "No sanity errors in ${INPUT_FILE} file" >> ${SANITY_RPT}
	echo "" >> ${SANITY_RPT}
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
# run sanity checks
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Running Preprocessor ${PREPROCESSOR}" >> ${LOG_DIAG}
if [ "${INPUT_FILE}" != "None" ]
then
    runSanityChecks
    STAT=$?
    checkStatus ${STAT} "runSanityChecks"
fi

${PREPROCESSOR}
STAT=$?
checkStatus ${STAT} "${PREPROCESSOR}"

#
# Run the load
#
echo "" >> ${LOG_DIAG}
date >> ${LOG_DIAG}
echo "Running loader ${LOADER}" >> ${LOG_DIAG}
${PYTHON} ${LOADER}
STAT=$?
checkStatus ${STAT} "${LOADER}"


#
# Do BCP
#
TABLE=MRK_Cluster

if [ -s "${OUTPUTDIR}/${TABLE}.bcp" ]
then

    echo "" >> ${LOG_DIAG}
    date >> ${LOG_DIAG}
    echo 'BCP data into MRK_Cluster'  >> ${LOG_DIAG}

    # BCP new data 
    ${PG_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}
fi

TABLE=MRK_ClusterMember

if [ -s "${OUTPUTDIR}/${TABLE}.bcp" ]
then

    echo "" >> ${LOG_DIAG}
    date >> ${LOG_DIAG}
    echo 'BCP data into MRK_ClusterMember'  >> ${LOG_DIAG}

    # BCP new data 
    ${PG_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}
fi

TABLE=ACC_Accession

if [ -s "${OUTPUTDIR}/${TABLE}.bcp" ]
then

    echo "" >> ${LOG_DIAG}
    date >> ${LOG_DIAG}
    echo 'BCP data into ACC_Accession'  >> ${LOG_DIAG}

    # BCP new data
    ${PG_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}
fi

TABLE=MGI_Property

if [ -s "${OUTPUTDIR}/${TABLE}.bcp" ]
then

    echo "" >> ${LOG_DIAG}
    date >> ${LOG_DIAG}
    echo 'BCP data into MGI_Property'  >> ${LOG_DIAG}

    # BCP new data
    ${PG_DBUTILS}/bin/bcpin.csh ${MGD_DBSERVER} ${MGD_DBNAME} ${TABLE} ${OUTPUTDIR} ${TABLE}.bcp ${COLDELIM} ${LINEDELIM} >> ${LOG_DIAG}
fi


cat - <<EOSQL | ${PG_DBUTILS}/bin/doisql.csh $0 >> ${LOG_DIAG}
select setval('mrk_cluster_seq', (select max(_Cluster_key) from MRK_Cluster));
select setval('mrk_clustermember_seq', (select max(_ClusterMember_key) from MRK_ClusterMember));
EOSQL

#
# run postload cleanup and email logs
#
shutDown
