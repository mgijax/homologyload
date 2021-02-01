#!/bin/sh

#
#
# invokes the Alliance Clustered preprocessor
#
#
# Usage:
#
#     preprocessAllianceClustered.sh
#

echo 'Running Preprocessor' >> ${LOG_DIAG}
${PYTHON} ${HOMOLOGYLOAD}/bin/preprocessAllianceClustered.py
STAT=$?
exit ${STAT}

exit 0
