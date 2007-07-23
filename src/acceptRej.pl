#!/usr/bin/perl -w

# acceptRej.pl
#
# Copyright (C) 2006   Naoki Takebayashi and Michael Hickerson
#
# This program is distributed with msBayes.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA


my $pdfOut = "figs.pdf";  # default pdf output file name

my $usage="Usage: $0 [-hn] [-p outPDF] [-s summary_stats] \n" .
          "                       [-t tolerance] obsData simData \n".
    "  -h: help\n".
    "  -n: print out the names of all available summary stats\n".
    "  -p: output pdf filename (default: $pdfOut)\n".
    "  -r: Simple rejection method without regression will be used\n".
    "  -s: statString (e.g. -s 'pi,wattTheta,pi.net,tajD.denom', <=default)\n".
    "      The summary statistics listed here will be used\n".
    "  -t: tolerance (a value between 0 an 1, default set in acceptRej.r)"
    ;


use Getopt::Std;
getopts('hndp:rt:s:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

our($opt_h, $opt_n, $opt_d, $opt_p, $opt_r, $opt_t, $opt_s);


use File::Copy;
use IO::File;
use POSIX qw(tmpnam);

# The names (not paths) of 3 R scripts required for this program.
my $mainRscript = "acceptRej.r";
my $make_pdRscript = "make_pd2005.r";
my $loc2plotRscript = "loc2plot.r";

# Adding the following paths to @INC, so we can find the R scripts.
# The R scripts should be in the same directory as this perl script,
# or inside of ../lib/msbayes/ relative to this perl script.
# i.e. If this script is inside of ~/bin/, the 3 r-scripts should be inside of
# ~/lib/msbayes/.
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib/msbayes";

# open a temporary file to store the dynamically created R script
do {$tmpR = tmpnam()} until $tmpRfh = 
    IO::File->new($tmpR, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpR) && -e $tmpR) {
	unlink($tmpR) || die "Couldn't unlink $tmpR : $!"
	}
};

# open a temp file to preprocess the observed data
do {$tmpObs = tmpnam()} until $tmpObsfh = 
    IO::File->new($tmpObs, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpObs) && -e $tmpObs) {
	unlink($tmpObs) || die "Couldn't unlink $tmpObs : $!"
	}
};


my ($simDat, $obsDat);
if(defined($opt_n)) {
    MkPrintNameRScript($tmpRfh);
} else {
    if (@ARGV != 2) {
	warn "ERROR: This script requires two arguments";
	die $usage;
    }
    
    ($obsDat, $simDat) = @ARGV;
    
    if (! (-r $simDat && -T $simDat )) {
	warn "ERROR: Problem with $simDat. Please give a readable, non-empty ".
	    "file name\n";
	die $usage;
    }
    if (! (-r $obsDat && -T $obsDat )) {
	warn "ERROR: Problem with $obsDat. Please give a readable, non-empty ".
	    "file name\n";
	die $usage;
    }

    if(defined($opt_p)) {
	$pdfOut=$opt_p;
    }
    CheckNBackupFile($pdfOut);
    MkStdAnalysisRScript($tmpRfh);
}

close($tmpRfh);  # the tmp R script is ready to use

if(defined($opt_d)) {
    open RSCRIPT, "<$tmpR" || die ("ERROR: Can't open $tmpR\n");
    warn "### debug: contents of R script file\n";
    while(<RSCRIPT>) {
	print STDERR $_;
    }
    close RSCRIPT;
    warn "### debug: end of R script file\n";
}

## preprocess the obs data
## Sometime the obs data isn't in 1 line, so converting this to 1 line (same
## format as the simulated data).
if (! defined($opt_n)) {
    open OBS, "<$obsDat" || die "Can't open $obsDat\n";
    while(<OBS>) {
	s/^\s+//; s/\s+$//;
	if ($. == 1) {
	    print $tmpObsfh $_;
	} else {
	    print $tmpObsfh "\t$_";
	}
    }
    close OBS;
    print $tmpObsfh "\n";
}

close $tmpObsfh;
# done with preproscess the obs data

# run R
my @output = `R --quiet --no-save --no-restore --slave < $tmpR`;

print @output;

exit (0);


# Making R script to print out stat names
sub MkPrintNameRScript {
    my $fh = shift;
    my $mainR = FindFile($mainRscript);
    if ($mainR eq '-1') {
	die "Can't find $mainRscript in directories:\n", join(":", @INC), "\n";
    }
    print $fh "source(\"$mainR\")\n";
    print $fh "printStatNames()\n";
    return;
}

# making the R script
sub MkStdAnalysisRScript {
    my $fh = shift;
    # find all required R scripts
    my $mainR = FindFile($mainRscript);
    if ($mainR eq '-1') {
	die "Can't find $mainRscript in directories:\n", join(":", @INC), "\n";
    }

    my $make_pd = FindFile($make_pdRscript);
    if ($make_pd eq '-1') {
	die "Can't find $make_pdRscript in directories:\n",join(":", @INC),"\n";
    }

    my $loc2plot = FindFile($loc2plotRscript);
    if ($make_pd eq '-1') {
	die "Can't find $loc2plotRscript in directories:\n",join(":", @INC),"\n";
    }

    print $fh "source(\"$make_pd\")\n";
    print $fh "source(\"$loc2plot\")\n";

    open MAIN_R_TMPL, "<$mainR";
    while(<MAIN_R_TMPL>) {
	s/^\s*source\s*\(\s*["']$make_pdRscript['"]\s*\)\s*$//;
	s/^\s*source\s*\(\s*["']$loc2plotRscript['"]\s*\)\s*$//;
	print $fh $_;
    }
    close MAIN_R_TMPL;
    
    # print $fh "source(\"$mainRt\")\n";
    print $fh "res <- stdAnalysis(\"$tmpObs\", \"$simDat\", pdf.outfile=\"$pdfOut\"";

    if (defined($opt_t)) {
	print $fh ", tol=$opt_t";
    }

    if (defined($opt_s)) {
	my $statString = MkStatString($opt_s);
	print $fh ", used.stats=c($statString)";
    }

    if (defined($opt_r)) {
	print $fh ", rejmethod=T";  # no regression
    } else {
	print $fh ", rejmethod=F";  # regression
    }
    print $tmpRfh ")\n";

    return;
}

sub MkStatString {
    my $ss = shift;
    $ss =~ s/^\s+//; $ss =~ s/\s+$//;
    if ($ss =~ /"/ || $ss =~ /'/) {  # ' or " is already included, do nothing
	return $ss;
    }
    $ss =~ s/\s*,\s*/","/g;
    $ss = "\"$ss\"";
    return $ss;
}

# This fucntion check if the argument (fileName) exists. If it exists,
# it get renamed to fileName.oldN, where N is a digit.
# In this way, no files will be overwritten.
sub CheckNBackupFile {
    my $fileName = shift;

    if (-e $fileName) {
        my $i = 1;
        while (-e "$fileName.old$i") {  # checking if the file exists
            $i++;
        }
        move("$fileName", "$fileName.old$i") ||
            die "Can't rename $fileName to $fileName.old$i";
    }
    # create the empty outfile, so other processes don't use the name.
    open(OUT,">$fileName");
    close(OUT);
}


# Check @INC library locations, and find a file
sub FindFile {
    my $name = shift;

    foreach my $dir (@INC) {
	if ( -e "$dir/$name" ) {
	    return "$dir/$name";
	}
    }
    return -1;
}
