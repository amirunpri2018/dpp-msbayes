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

ARCHITECTURES="m32 m64"
 
# number of cpus to use during compile
COMPILETHREADS=2
 
COMPILER=gcc
 
# configure flags
CONFFLAGS="--disable-shared"
 
for ARCH in $ARCHITECTURES
do
    echo Compiling for ${ARCH}...
     
    # copy source to new directory
    ARCH_DIR="${GSLVER}_${ARCH}"
    INSTALL_DIR="${BUILD_DIR}/${ARCH}"
    cp -r "${GSLVER}" "$ARCH_DIR"
    cd "$ARCH_DIR"
     
    env CC="$COMPILER" CFLAGS="-$ARCH" LDFLAGS="-$ARCH" ./configure $CONFFLAGS --prefix="$INSTALL_DIR"
     
    # compile for architecture
    make clean
    nice make -j $COMPILETHREADS
    make install
     
    echo
    echo
    echo ${ARCH} headers and binaries are in:
    echo "    $INSTALL_DIR"

    env_path="${BUILD_DIR}/gsl-${ARCH}-env.sh"
    echo "#!/bin/sh" > "$env_path"
    echo export PATH="${BUILD_DIR}/${ARCH}/bin:${PATH}" >> "$env_path"
    echo export LD_LIBRARY_PATH="${BUILD_DIR}/${ARCH}/lib:${LD_LIBRARY_PATH}" >> "$env_path"
    echo export PKG_CONFIG_PATH="${BUILD_DIR}/${ARCH}/lib/pkgconfig:${PKG_CONFIG_PATH}" >> "$env_path"
    echo export GSL_PREFIX="${BUILD_DIR}/${ARCH}" >> "$env_path"

    cd "$BASE_DIR"

done
