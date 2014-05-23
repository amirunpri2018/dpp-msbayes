#!/usr/bin/env bash
 
# Modified from script posted by Mark Dayel to
# <http://www.dayel.com/blog/2009/09/05/building-gsl-universal-binary/> on
# 2009-09-05.

GSL_TAR_BALL="gsl-1.16.tar.gz"
BASE_DIR="$(pwd)"
BUILD_DIR="${BASE_DIR}/build-gsl"
if [ -d "$BUILD_DIR" ]
then 
    echo "ERROR: build directory '$BUILD_DIR' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "$BUILD_DIR"
fi
SOURCEFILE="${BUILD_DIR}/${GSL_TAR_BALL}"
curl -o "$SOURCEFILE" "http://mirror.sdunix.com/gnu/gsl/${GSL_TAR_BALL}"
GSLVER=${SOURCEFILE%.tar.gz}
tar -xzf "$SOURCEFILE" -C "$BUILD_DIR"

# possible archs: i386 ppc x86_64 ppc64
# note: ppc64 doesn't work in Snow Leopard
ARCHITECTURES="i386 ppc x86_64"
 
# number of cpus to use during compile
COMPILETHREADS=2
 
COMPILER=gcc
#COMPILER=gcc-4.2
#COMPILER=/Developer/usr/bin/llvm-g++-4.2
 
# configure flags
CONFFLAGS="--disable-shared"
 
# compiler flags
COMPFLAGS="-O3"
 
 
for ARCH in $ARCHITECTURES
do
    echo Compiling for ${ARCH}...
     
    # copy source to new directory
    cp -r "${GSLVER}" "${GSLVER}_$ARCH"
    cd "${GSLVER}_$ARCH"
     
    # configure for architecture
    if [ "$ARCH" != "ppc64" ] ; then  # not sure why we have to do cross compile for ppc64 but not ppc...
        env CC="$COMPILER" CFLAGS="-arch $ARCH $COMPFLAGS" LDFLAGS="-arch $ARCH" ./configure $CONFFLAGS --build=$ARCH
    else
        env CC="$COMPILER" CFLAGS="-arch $ARCH $COMPFLAGS" LDFLAGS="-arch $ARCH" ./configure $CONFFLAGS --host=$ARCH
    fi
     
    # compile for architecture
    make clean ; nice make -j $COMPILETHREADS
     
    cd "$BASE_DIR"
     
done
 
INCLUDE_DIR="${BUILD_DIR}/include"
LIB_DIR="${BUILD_DIR}/libs"
mkdir -p "$INCLUDE_DIR"
mkdir -p "$LIB_DIR"
cp `find ${GSLVER} -name "*.h"` "$INCLUDE_DIR"
 
echo
echo
echo Making Universal Binaries...
 
find ${GSLVER}_$ARCH -name "*.a"
 
for LIBRARY in `find ${GSLVER}_$ARCH -name "*.a"`
do
    echo
    LIB=`basename $LIBRARY`
    echo Library: $LIB
    find . -name $LIB
    lipo -create `find . -name $LIB | xargs` -output "${LIB_DIR}/${LIB}"
    lipo -info "${LIB_DIR}/${LIB}"
done
 
echo
echo Universal Binaries:
echo
lipo -info "$LIB_DIR"/*.a

