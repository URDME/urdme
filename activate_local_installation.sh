#!/bin/bash

####################################################
DEFAULT=`pwd`"/urdme"
####################################################
if [ ! -d $DEFAULT ]; then
    echo "Error: $DEFAULT not found, exiting"
    exit
fi
####################################################
echo "This script will activated $DEFAULT as the current URDME installation on this computer."
echo "Press return to continue,  or ^C (ctrl + c) to EXIT"
read -s
####################################################
URDME_ROOT=$DEFAULT
####################################################
HAS_MATLAB=`which matlab`
if [ "$HAS_MATLAB" = "" ]; then
    echo "ERROR: Matlab executable not found in your PATH"
    exit
fi
####################################################
echo -n "URDME version: "
VER=`./urdme/bin/urdme_init -v`
echo $VER
###
echo -n "Matlab: "
echo `./urdme/bin/urdme_init -m`
###
echo -n "Matlab Arch: "
echo `./urdme/bin/urdme_init -a`
####################################################
echo "URDME installation destination: $URDME_ROOT"
echo "To specify a different use the 'install.sh' script'"
####################################################
####################################################
####################################################
echo "Press return to continue with installation, or ^C (ctrl + c) to EXIT"
read -s
####################################################
if [ ! -d /usr/local/bin ]; then
    mkdir -p /usr/local/bin
fi
rm -f /usr/local/bin/urdme_init
ln -s "$URDME_ROOT/bin/urdme_init" /usr/local/bin
####################################################
echo "Press return to have the install scripts add the URDME folders to the Matlab path, or ^C (ctrl + c) to EXIT"
read -s
####################################################
./set_matlab_path.pl $URDME_ROOT >  /dev/null
####################################################
echo "Local URDME installation activated"
####################################################


