#!/bin/sh

#
#
# invokes the Alliance direct preprocessor
#
#
# Usage:
#
#     preprocessAllianceDirect.sh
#

echo 'Running Preprocessor' >> ${LOG_DIAG}
${PYTHON} ${HOMOLOGYLOAD}/bin/preprocessAllianceDirect.py
STAT=$?
exit ${STAT}

exit 0
