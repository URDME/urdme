#!/bin/bash

####################################################
DEFAULT="/usr/local/urdme"
if [ "$1" == "" ]; then
    URDME_ROOT=$DEFAULT
else
    URDME_ROOT=$1
fi
####################################################
# if $1 is "PREFIX=/blah/blah/"
####################################################
HAS_MATLAB=`which matlab`
if [ "$HAS_MATLAB" = "" ]; then
    echo "ERROR: Matlab executable not found in your PATH"
    exit
fi
####################################################
echo -n "URDME version: "
echo `./bin/urdme_init -v`
###
echo -n "Matlab: "
echo `./bin/urdme_init -m`
###
echo -n "Matlab Arch: "
echo `./bin/urdme_init -a`
####################################################
echo "URDME installation destination: $URDME_ROOT"
echo "To specify a different installation location, exit and call this script with the location as the argument."
echo "    install.sh /path/to/install/directory/ "
if [ ! -d "$URDME_ROOT" ]; then
    mkdir -p $URDME_ROOT
else
    FILE_COUNT=`/bin/ls -l $URDME_ROOT | wc -l`
    if [ "$FILE_COUNT" != "0" ]; then
        echo "$URDME_ROOT already exists and is not empty.  Please remove all previous installations of URDME."
        exit
    fi
    echo "FILE_COUNT=$FILE_COUNT"
fi
if [ ! -w "$URDME_ROOT" ]; then
    echo "ERROR: '$URDME_ROOT' is not writeable"
    echo "Run this script with the correct permissions."
    exit
fi
####################################################
####################################################
####################################################
echo "Press return to continue with installation, or ^C (ctrl + c) to EXIT"
read -s
####################################################
cd ..
cp -r ./urdme/* "$URDME_ROOT"
cd urdme
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
echo "URDME installation complete"
####################################################
