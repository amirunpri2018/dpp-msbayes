Table of Contents
=================

 -  [Installation](#installation)
    -  [Requirements](#requirements)
    -  [Quick Default Install](#quick-default-install)
    -  [More Install Details](#more-install-details)
        -  [Specifying GSL location](#specifying-gsl-location)
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

Quick Default Install
---------------------

Clone the code repository and `cd` into it:

    $ git clone https://github.com/joaks1/dpp-msbayes.git
    $ cd dpp-msbayes

Initiate and clone the `ABACUS` submodule:

    $ git submodule init
    $ git submodule update

I like to create a new directory for compiling, so that if all goes wrong (or
after all goes well) I can blow it away:

    $ mkdir build
    $ cd build

Call the `configure` script to generate the `Makefile`, and compile with
`make`:

    $ ../configure
    $ make

If you have the `Check` package installed and run the test suite:

    $ make check

Finally, install to the default location (`/usr/local/bin` and
`/usr/local/lib`; you will most likely need to use `sudo` to install to
`/usr/local`):

    $ sudo make install

More Install Details
--------------------

All the code snippets below assume you created and moved into a scratch
directory created at the top level of the `dpp-msbayes` package, like so:

    $ cd dpp-msbayes
    $ mkdir build
    $ cd build

You can peruse all configuration, build, and install options via:

    $ ../configure -h

Let's take a look at some of the more important `configure` options.

### Specifying GSL location

Part of the output of `../configure -h` is:

    --with-gsl=yes/no/DIR   if yes, configure will look for a system-level GSL
                            installation (this is the default); if no, configure
                            will only build and install gsl-independent ABACUS
                            tools; if a directory is provided, configure will
                            use the gsl library installed at that location. For
                            example, if the non-system-level GSL library files
                            you wish to use are in /home/joebob/lib, then use
                            --with-gsl=/home/joebob.

If GSL is installed in a rational way, the configure script should be able to
locate it by default via `gsl-config` or `pkg-config`. However, if you have
a non-standard install of GSL, and the configure script cannot find it by
default, simply tell it where it is:

    $ ../configure --with-gsl=/path/to/where/gsl/lives

### Changing the install location

If you prefer not to install to the default `/usr/local/` location,
you can easily change this using the `--prefix` flag:

    $ ../configure --prefix=${HOME}/Environment

On my laptop, this would install the binaries and executables to
`/home/jamie/Environment/bin` and some R libraries to
`/home/jamie/Environment/lib`. If you do not have admin privileges on your
machine you can use this option to install the software anywhere within your
home folder. If you do specify a location within your home folder you can
install without using "super user" privileges:

    $ make
    $ make check  # if you have Check installed
    $ make install

Note the lack of `sudo` in `make install` that was used above.

