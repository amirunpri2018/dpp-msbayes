Table of Contents
=================

 -  [Installation](#installation)
    -  [Requirements](#requirements)
    -  [Quick Default Install](#quick-default-install)
    -  [More Install Details](#more-install-details)
        -  [Changing the install location](#changing-the-install-location)

Installation
============

I strongly recommend you use this software via the Python multiprocessing
package `PyMsBayes` (<https://github.com/joaks1/PyMsBayes>). This package comes
bundled with pre-compiled `dpp-msbayes` executables (and the original `msBayes`
executable, as well, which you can use interchangeably). By using `PyMsBayes`,
you can avoid compiling from source code and using the tools in this package
"by hand."

Requirements
------------

You will need a C compiler and the GNU Scientific Library (GSL;
<http://www.gnu.org/software/gsl/>). If you want to run the test suite, you
will also need the `Check` package (<http://check.sourceforge.net/>).

This version of `msBayes` now uses CMake (<http://www.cmake.org/>) for the
build system. If your computer does not have `cmake`, you will need to install
it from <http://www.cmake.org/cmake/resources/software.html>. You can test if
you have `cmake` by running the following on the command line:

    $ cmake --version

Quick Default Install
---------------------

Clone the code repository and `cd` into it:

    $ git clone https://github.com/joaks1/dpp-msbayes.git
    $ cd dpp-msbayes

To build the binaries:

    $ ./build.sh

Finally, install to the default location (`/usr/local/bin` and
`/usr/local/lib`; you will most likely need to use `sudo` to install to
`/usr/local`):

    $ cd build
    $ sudo make install

You can also build and install in one step if you are impatient. Instead of the
three previous commands, do:

    $ sudo ./build.sh --install

More Install Details
--------------------

All the code snippets below assume you are in the top level of the
`dpp-msbayes` package.

You can peruse all configuration, build, and install options via:

    $ ./build.sh -h

Let's take a look at some of the more important `build.sh` options.

### Changing the install location

If you prefer not to install to the default `/usr/local/` location, you can
easily change this using the `-p|--prefix` flag:

    $ ./build.sh --prefix ${HOME}/Environment --install

On my laptop, this would install the binaries and executables to
`/home/jamie/Environment/bin` and some R libraries to
`/home/jamie/Environment/lib`. If you do not have admin privileges on your
machine you can use this option to install the software anywhere within your
home folder. If you do specify a location within your home folder you can
install without using "super user" privileges (`sudo`).

