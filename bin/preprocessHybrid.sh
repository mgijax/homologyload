#!/bin/sh

#
#
# invokes the Hybrid preprocessor
#
#
# Usage:
#
#     preprocessHybrid.sh
#

echo 'Running Preprocessor' >> ${LOG_DIAG}
${PYTHON} ${HOMOLOGYLOAD}/bin/preprocessHybrid.py
STAT=$?
exit ${STAT}

exit 0
