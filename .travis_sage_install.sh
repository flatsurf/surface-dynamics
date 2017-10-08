#!/usr/bin/env bash
# Sage installation script for Travis-CI

cd ${HOME}
wget ${SAGE_MIRROR}/${SAGE_BINARY}
tar xf ${SAGE_BINARY}
mv SageMath ${SAGE_ROOT}
ls ${SAGE_ROOT}

if [ ! -x "${SAGE}" ] ;
then
    echo "Installation failed!"
    exit 1
fi

# the first time we launch Sage the relocalization script is run
${SAGE} -c ''
if [ $? -ne 0 ]
then
    echo "Sage does not start"
    exit 1
fi
