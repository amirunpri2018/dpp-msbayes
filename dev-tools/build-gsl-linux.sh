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
ARCHITECTURES="m32 m64"
 
# number of cpus to use during compile
COMPILETHREADS=2
 
COMPILER=gcc
#COMPILER=gcc-4.2
#COMPILER=/Developer/usr/bin/llvm-g++-4.2
 
# configure flags
CONFFLAGS="--disable-shared"
 
# compiler flags
# COMPFLAGS="-O3"
 
 
for ARCH in $ARCHITECTURES
do
    echo Compiling for ${ARCH}...
     
    # copy source to new directory
    ARCH_DIR="${GSLVER}_${ARCH}"
    cp -r "${GSLVER}" "$ARCH_DIR"
    cd "$ARCH_DIR"
     
    # configure for architecture
    # if [ "$ARCH" != "ppc64" ] ; then  # not sure why we have to do cross compile for ppc64 but not ppc...
    #     env CC="$COMPILER" CFLAGS="-arch $ARCH $COMPFLAGS" LDFLAGS="-arch $ARCH" ./configure $CONFFLAGS --build=$ARCH
    # else
    env CC="$COMPILER" CFLAGS="-$ARCH" LDFLAGS="-$ARCH" ./configure $CONFFLAGS
    # fi
     
    # compile for architecture
    make clean ; nice make -j $COMPILETHREADS
     
    echo
    echo
    echo Copying ${ARCH} headers and binaries...
    INCLUDE_DIR="${BUILD_DIR}/${ARCH}/include"
    LIB_DIR="${BUILD_DIR}/${ARCH}/libs"
    mkdir -p "$INCLUDE_DIR"
    mkdir -p "$LIB_DIR"
    cp `find ${ARCH_DIR} -name "*.h"` "$INCLUDE_DIR"
    cp `find ${ARCH_DIR} -name "*.a"` "$LIB_DIR"
     
    echo
    echo
    echo ${ARCH} headers are in:
    echo "    $INCLUDE_DIR"
    
    echo
    echo
    echo ${ARCH} static libraries:
    file "${LIB_DIR}"/*.a

    cd "$BASE_DIR"
     
done
