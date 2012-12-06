#
#  Creates the core URDME distribution. 
#
#
# A. Hellander 2010-03-21
# B. Drawert   2010-06-01
VER=`urdme_init -v`
DEST="../urdme-$VER"
echo "Creating distribution $DEST"

mkdir $DEST
mkdir $DEST/urdme 
mkdir $DEST/urdme/bin
mkdir $DEST/urdme/build
mkdir $DEST/urdme/src
mkdir $DEST/urdme/src/nsm
mkdir $DEST/urdme/msrc
mkdir $DEST/urdme/msrc/fvm
mkdir $DEST/urdme/include
mkdir $DEST/urdme/comsol
mkdir $DEST/examples/
mkdir $DEST/examples/mincde
mkdir $DEST/examples/bistab
mkdir $DEST/urdme/doc

# distribution files
cp distribute/AUTHORS $DEST
cp distribute/README $DEST
cp distribute/COPYING $DEST
cp distribute/install.sh $DEST
cp distribute/set_matlab_path.pl $DEST

# Makefiles
cp build/Makefile.nsm $DEST/urdme/build/
cp bin/urdme_init $DEST/urdme/bin/

# NSM Solver files
cp  src/nsm/nsm.h $DEST/urdme/src/nsm
cp  src/nsm/nsm.c $DEST/urdme/src/nsm
cp  src/nsm/nsmcore.c $DEST/urdme/src/nsm
cp  src/nsm/binheap.c $DEST/urdme/src/nsm
cp  src/nsm/binheap.h $DEST/urdme/src/nsm

# General source files
cp  src/report.c   $DEST/urdme/src
cp  src/matmodel.c $DEST/urdme/src
cp  src/read_matfile.c $DEST/urdme/src

# Matlab interface
cp  msrc/urdme.m  $DEST/urdme/msrc
cp  msrc/urdme_compile.m  $DEST/urdme/msrc
cp  msrc/rdme2mat.m  $DEST/urdme/msrc
cp  msrc/urdme_startup.m $DEST/urdme/msrc
cp  msrc/urdme_validate.m $DEST/urdme/msrc
cp  msrc/urdme_addsol.m $DEST/urdme/msrc
cp  msrc/urdme_savesol.m $DEST/urdme/msrc

# headers
cp  include/matmodel.h     $DEST/urdme/include
cp  include/report.h       $DEST/urdme/include
cp  include/propensities.h $DEST/urdme/include
cp  include/read_matfile.h $DEST/urdme/include

# comsol interface routines
cp  comsol/comsol2urdme.m $DEST/urdme/comsol
cp  comsol/urdme2comsol.m $DEST/urdme/comsol

# mincde example
cp  examples/mincde/mincde.m $DEST/examples/mincde
cp  examples/mincde/mincde.c $DEST/examples/mincde
cp  examples/mincde/coli.mph $DEST/examples/mincde
cp  examples/mincde/temporal_average.m $DEST/examples/mincde

# documentation
cp  doc/manual.pdf $DEST/urdme/doc

# create tarball
#tar -cf urdme.tar urdme
#gzip urdme.tar
#rm -rf urdme



