#!/bin/sh

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-p|--prefix <INSTALL-PREFIX>] [-s|--static] [-i|--install] [-t|--test]"
    echo "  -h|--help       show help message and exit."
    echo "  -s|--static     build statically linked binaries."
    echo "  -p|--prefix     install path. Default: '/usr/local'."
    echo "  -t|--test       run the test suite after building."
    echo "  -i|--install    install the executables."
    echo "  --old           include old model in build."
    echo ""
}

build_dpp_msbayes () {
    cmake_args=$@
    if [ -n "$DPP_MSBAYES_BASE_DIR" ]
    then
        base_dir="$DPP_MSBAYES_BASE_DIR"
    else
        echo "ERROR: DPP_MSBAYES_BASE_DIR was not defined prior to calling "
        echo "build_dpp_msbayes."
        exit 1
    fi
    if [ -n "$DPP_MSBAYES_BUILD_DIR" ]
    then
        build_dir="$DPP_MSBAYES_BUILD_DIR"
    else
        build_dir="${base_dir}/build"
    fi
    if [ -d "$build_dir" ]
    then
        echo "ERROR: build directory '$build_dir' already exists."
        echo "To reconfigure, please remove this directory and re-run this script."
        exit 1
    else
        mkdir -p "$build_dir"
    fi
    
    # configure make files and build
    cd "$build_dir"
    echo "Configuring make files..."
    echo "    "${base_dir}"/ $cmake_args | xargs cmake"
    echo "${base_dir}/ $cmake_args" | xargs cmake || exit 1
    echo "Building..."
    make || exit 1
    if [ -n "$DPP_MSBAYES_RUN_TESTS" ]
    then
        echo "Building and running test suite..."
        make check || exit 1
    fi
    if [ -n "$DPP_MSBAYES_RUN_INSTALL" ]
    then
        echo "Installing..."
        make install || exit 1
    fi
    cd "$base_dir"
    echo "Done!"
}


# get location of script
DPP_MSBAYES_BASE_DIR=""
this_dir=`dirname "$0"`
if [ "$this_dir" = "." ]
then
    export DPP_MSBAYES_BASE_DIR="$(pwd)"
else
    cd "$this_dir"
    export DPP_MSBAYES_BASE_DIR="$(pwd)"
fi

# make sure abacus is here and up to date
git submodule init
git submodule update

# process args
extra_args=""
static=""
install_prefix=""
DPP_MSBAYES_RUN_TESTS=""
DPP_MSBAYES_RUN_INSTALL=""
include_old=""
universal_mac_build=""
universal_linux_build=""
while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -p| --prefix)
            shift
            install_prefix=$1
            ;;
        -s| --static)
            static=1
            ;;
        -t| --test)
            DPP_MSBAYES_RUN_TESTS=1
            export DPP_MSBAYES_RUN_TESTS
            ;;
        -i| --install)
            DPP_MSBAYES_RUN_INSTALL=1
            export DPP_MSBAYES_RUN_INSTALL
            ;;
        --old)
            include_old=1
            ;;
        --universal-mac-build)
            universal_mac_build=1
            include_old=1
            static=1
            ;;
        --universal-linux-build)
            universal_linux_build=1
            include_old=1
            static=1
            ;;
        * )
            extra_args="$extra_args $1"
    esac
    shift
done
args=""
if [ -n "$static" ]
then
    args="${args} -DSTATIC_LINKING=YES"
fi
if [ -n "$include_old" ]
then
    args="${args} -DINCLUDE_OLD_MODEL=YES"
fi
if [ -n "$extra_args" ]
then
    args="${args} ${extra_args}"
fi
if [ -n "$universal_linux_build" ]
then
    for arch in "32" "64"
    do
        args="${args} -DCMAKE_INSTALL_PREFIX=${DPP_MSBAYES_BASE_DIR}/install/m${arch}"
        args="${args} -DCMAKE_C_FLAGS=-m${arch} -DCMAKE_LD_FLAGS=-m${arch}"
        export DPP_MSBAYES_BUILD_DIR="${DPP_MSBAYES_BASE_DIR}/build/m${arch}"
        build_dpp_msbayes $args
    done
fi
if [ -n "$install_prefix" ]
then
    args="${args} -DCMAKE_INSTALL_PREFIX=${install_prefix}"
fi
if [ -n "$universal_mac_build" ]
then
    args="${args} -DCMAKE_C_FLAGS=\"-arch ppc -arch i386 -arch x86_64\""
fi

# check for build directory
export DPP_MSBAYES_BUILD_DIR="${DPP_MSBAYES_BASE_DIR}/build"
build_dpp_msbayes $args
