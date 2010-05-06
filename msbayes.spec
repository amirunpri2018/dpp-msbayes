%define debug_package  %{nil}

%define name msbayes
%define version 20100505
%define release 1

%define prefix   /usr/local
%define builddir $RPM_BUILD_DIR/%{name}-%{version}

Summary: A program for testing comparative phylogeographic histories
Name: %{name}
Version: %{version}
Release: %{release}
Group: Applications/Scientific
License: GPL
Packager: Naoki Takebayashi <ffnt@uaf.edu>
URL: http://msbayes.sourceforge.net/
Source0: http://msbayes.sourceforge.net/msbayes/%{name}-%{version}.tgz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root
Requires: R
BuildRequires: gsl-devel gsl-static

%description

msBayes is a pipeline for testing comparative phylogeographic
histories using hierarchical ABC, developped by Michael Hickerson.  It
can be used to test for simultaneous divergence (TSD) of multiple
taxon pairs.

For the full description, see:

2006 Hickerson, M.J., E, Stahl, and H.A. Lessios. Test for
simultaneous divergence (TSD) using approximate Bayesian computation
(ABC). In press (Evolution)

The binaries are statically linked, so you do not need to install GSL
if you use this RPM.

%prep
%setup -n %{name}-%{version}

%build
cd src; make LDFLAGS=-static

%install
function CheckBuildRoot() {
    # do a few sanity checks on the BuildRoot
    # to make sure we don't damage a system
    case "${RPM_BUILD_ROOT}" in
        ''|' '|/|/bin|/boot|/dev|/etc|/home|/lib|/mnt|/root|/sbin|/tmp|/usr|/var)
            echo "Yikes!  Don't use '${RPM_BUILD_ROOT}' for a BuildRoot!"
            echo "The BuildRoot gets deleted when this package is rebuilt;"
            echo "something like '/tmp/build-blah' is a better choice."
            return 1
            ;;
        *)  return 0
            ;;
    esac
}

function CleanBuildRoot() {
    if CheckBuildRoot; then
	rm -rf "${RPM_BUILD_ROOT}"
    else
        exit 1
    fi
}

CleanBuildRoot

cd src; make PREFIX=${RPM_BUILD_ROOT}/%{prefix} install

%clean
rm -r $RPM_BUILD_ROOT
rm -r %{builddir}

%files
%defattr(-,root,root)
%doc documents/* INSTALL LICENSE README VERSION
%{prefix}/bin/*
%{prefix}/lib/msbayes

%changelog
* Wed May 05 2010 Naoki Takebayashi <ffnt@uaf.edu> [20100505-1]
- version update
- multi-locus version

* Thu Nov 06 2008 Naoki Takebayashi <ffnt@uaf.edu> [20081106-1]
- version update

* Fri Jun 13 2008 Naoki Takebayashi <ffnt@uaf.edu> [20080613-1]
- version update

* Thu May 15 2008 Naoki Takebayashi <ffnt@uaf.edu> [20080515-1]
- version update

* Mon Nov  6 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061106-1]
- version update

* Fri Oct 13 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061003-2]
- compiling with static library.

* Fri Oct 13 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061005-1]
- new version with correct example files

* Tue Oct 03 2006 Naoki Takebayashi <ffnt@uaf.edu> [20061003-1]
- first release
