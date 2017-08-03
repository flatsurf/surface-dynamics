#!/usr/bin/env bash
# Sage installation script for Travis-CI

if [ ! -x "${SAGE}" ] ;
then
    cd ${HOME}
    echo "Downloading and installing Sage..."
    wget ${SAGE_MIRROR}/${SAGE_BINARY}
    tar xf ${SAGE_BINARY}
    echo "done"

    if [ ! -x "${SAGE}" ] ;
    then
        echo "Installation failed!"
        exit 1
    fi

    # the first time we launch Sage the relocalization script is run
    ${SAGE} -c ''
fi
