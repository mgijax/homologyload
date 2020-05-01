#!/bin/sh

#
#
# invokes the HGNC preprocessor
#
#
# Usage:
#
#     preprocessHGNC.sh
#

echo 'Running Preprocessor' >> ${LOG_DIAG}
${PYTHON} ${HOMOLOGYLOAD}/bin/preprocessHGNC.py
STAT=$?
exit ${STAT}

exit 0
