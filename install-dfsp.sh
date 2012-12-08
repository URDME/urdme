#!/bin/bash

####################################################
DEFAULT="/usr/local/urdme"
if [ "$1" == "" ]; then
    if [ -e /usr/local/bin/urdme_init ]; then
        URDME_ROOT=`/usr/local/bin/urdme_init -r`
    else
        URDME_ROOT=$DEFAULT
    fi
else
    URDME_ROOT=$1
fi
echo "URDME installation destination: $URDME_ROOT"
echo "To specify a different installation location, exit and call this script with the location as the argument."
echo "    install.sh /path/to/install/directory/ "
####################################################
if [ ! -d "$URDME_ROOT" ]; then
    echo "ERROR: '$URDME_ROOT' not found"
    echo "ERROR: can not find core installation of URDME.  Please install the core distribution or confirm that the installtion directory is correct .";
fi
if [ ! -w "$URDME_ROOT" ]; then
    echo "ERROR: '$URDME_ROOT' is not writeable"
    echo "Run this script with the correct permissions."
    exit
fi
####################################################
if [ ! -e "$URDME_ROOT/bin/urdme_init" ]; then
    echo "ERROR: '$URDME_ROOT/bin/urdme_init' not found.";
    echo "ERROR: can not find core installation of URDME.  Please install the core distribution or confirm that the installtion directory is correct .";
    exit;
fi
####################################################
####################################################
####################################################
echo "Press return to continue with installation, or ^C (ctrl + c) to EXIT"
read -s
####################################################
cp -r ./urdme/* "$URDME_ROOT"
####################################################
echo "DFSP-URDME installation complete"
####################################################


