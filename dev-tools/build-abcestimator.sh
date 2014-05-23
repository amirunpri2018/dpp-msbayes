#!/usr/bin/env bash
 
ABC_ZIP="ABCtoolbox.zip"
BASE_DIR="$(pwd)"
BUILD_DIR="${BASE_DIR}/build-abcestimator"
if [ -d "$BUILD_DIR" ]
then 
    echo "ERROR: build directory '$BUILD_DIR' already exists."
    echo "To recompile, please remove this directory and re-run this script."
    exit 1
else
    mkdir -p "$BUILD_DIR"
fi
ABC_ZIP_PATH="${BUILD_DIR}/${ABC_ZIP}"
ABC_BASE_DIR="${ABC_ZIP_PATH%.zip}"
curl -o "$ABC_ZIP_PATH" "http://cmpg.unibe.ch/software/ABCtoolbox/${ABC_ZIP}"
unzip "$ABC_ZIP_PATH" -d "$BUILD_DIR"

cd "${ABC_BASE_DIR}/source/ABCestimator"

if [ "$(uname)" = "Darwin" ]
then
    MATH_LIB="$(find /usr -name "libm.dylib")"
    g++ -arch ppc -arch i386 -arch x86_64  -Wl,-search_paths_first -Wl,-headerpad_max_install_names -Bstatic *.cpp  -o ABCestimator "$MATH_LIB"
else
    MATH_LIB="$(find /usr -name "libm.a")"
    g++ -m32  -Wl,-search_paths_first -Wl,-headerpad_max_install_names -Bstatic *.cpp  -o ABCestimator32 "$MATH_LIB"
    g++ -m64  -Wl,-search_paths_first -Wl,-headerpad_max_install_names -Bstatic *.cpp  -o ABCestimator64 "$MATH_LIB"
fi

cd "$BASE_DIR"
