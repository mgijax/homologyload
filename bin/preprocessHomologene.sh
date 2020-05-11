#!/bin/sh

#
#
# invokes the Homologene preprocessor
#
#
# Usage:
#
#     preprocessHomologene.sh
#

echo 'Running Preprocessor' >> ${LOG_DIAG}
${PYTHON} ${HOMOLOGYLOAD}/bin/preprocessHomologene.py
STAT=$?
exit ${STAT}

exit 0
