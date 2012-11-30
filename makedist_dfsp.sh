#
#  Adds the DFSP solver to the core URDME distribution. 
#
#
# B. Drawert   2010-06-01
VER=`urdme_init -v`
DEST="../urdme-$VER"
echo "Creating distribution $DEST"

mkdir $DEST
mkdir $DEST/urdme 
mkdir $DEST/urdme/build
mkdir $DEST/urdme/src
mkdir $DEST/urdme/src/dfsp
mkdir $DEST/urdme/msrc
mkdir $DEST/examples/
mkdir $DEST/examples/polarization

# Makefiles
cp build/Makefile.dfsp $DEST/urdme/build/

# distribution files
if [ ! -e "$DEST/install.sh" ]; then
    echo "$DEST/install.sh not found, copying install-dfsp.sh"
    cp distribute/AUTHORS $DEST
    cp distribute/README $DEST
    cp distribute/COPYING $DEST
    cp distribute/install-dfsp.sh $DEST/install.sh
else
    echo "$DEST/install.sh found"
fi

# DFSP Solver files
cp  src/dfsp/dfsp.h $DEST/urdme/src/dfsp
cp  src/dfsp/dfsp.c $DEST/urdme/src/dfsp
cp  src/dfsp/dfspcore.c $DEST/urdme/src/dfsp
cp  src/dfsp/dfsp_diffusion.c $DEST/urdme/src/dfsp
cp  src/dfsp/dfsp_reactions.c $DEST/urdme/src/dfsp

# General source files
  
# Matlab interface
cp  msrc/urdme_init_dfsp.m $DEST/urdme/msrc

# headers

# comsol interface routines

# polarization example
cp  examples/polarization/polarization.m $DEST/examples/polarization
cp  examples/polarization/polarization.c $DEST/examples/polarization
cp  examples/polarization/polarization.mph $DEST/examples/polarization


# documentation

# create tarball
#tar -cf urdme.tar urdme
#gzip urdme.tar
#rm -rf urdme



