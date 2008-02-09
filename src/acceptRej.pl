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
my $defaultTolerance = 0.002;  # default tolerance
my $statString = 'pi,wattTheta,pi.net,tajD.denom';  # default stat strings

my $usage="Usage: $0 [-hnacr] [-p outPDF] [-s summary_stats] \n" .
          "                       [-t tolerance] obsData simData \n".
    "  -h: help\n".
    "  -n: print out the names of all available summary stats\n".
    "  -p: output pdf filename (default: $pdfOut)\n".
    "  -r: Simple rejection method without regression will be used\n".
    "  -s: statString (e.g. -s '$statString', <=default)\n".
    "      The summary statistics listed here will be used\n".
    "  -a: old analysis with all R processing (nobody needs it)\n" .
    "  -t: tolerance (a value between 0 an 1, default: $defaultTolerance)\n" .
    "  -c: analysis based on prior constrained by Psi=2 to allow obtaining\n".
    "      posteriors of extra hyper-parameters Psi1, Psi2, tau1, and tau2 \n"
    ;


use Getopt::Std;
getopts('ahncdp:rt:s:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));

our($opt_a, $opt_h, $opt_n, $opt_d, $opt_p, $opt_r, $opt_t, $opt_s, $opt_c);

use File::Copy;
use IO::File;
use POSIX qw(tmpnam);

# The names (not paths) of 3 R scripts and a C prog required for this program.
my $mainRscript = "acceptRej.r";
my $make_pdRscript = "make_pd2005.r";
my $loc2plotRscript = "loc2plot.r";
my $rejectionExe = "msReject"; # rejection program

# Adding the following paths to @INC, so we can find the R scripts.
# The R scripts should be in the same directory as this perl script,
# or inside of ../lib/msbayes/ relative to this perl script.
# i.e. If this script is inside of ~/bin/, the 3 r-scripts should be inside of
# ~/lib/msbayes/.
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib/msbayes";
use lib ".";  # looks for Required files in the current directory at first
              # -d will print out the search path
if (defined($opt_d)) {
    print STDERR "FILEINFO: searching path is ";
    print STDERR join(":", @INC);
    print STDERR "\n";
}
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

# open a temp file to preprocess the data with msrejection
my $tmpSimDatfh;
do {$tmpSimDat = tmpnam()} until $tmpSimDatfh = 
    IO::File->new($tmpSimDat, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpSimDat) && -e $tmpSimDat) {
	unlink($tmpSimDat) || die "Couldn't unlink $tmpSimDat : $!"
	}
};

# open a temp file to extract the prior columns.
do {$tmpPrior = tmpnam()} until $tmpPriorfh = 
    IO::File->new($tmpPrior, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpPrior) && -e $tmpPrior) {
	unlink($tmpPrior) || die "Couldn't unlink $tmpPrior : $!"
	}
};

if (defined($opt_s)) {
    $statString = MkStatString($opt_s);
}

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


{   # Making sure R is installed
    my $checkR = system("which R > /dev/null");
    die "\nERROR: Can't find R.  Make sure R is installed and in your PATH.\n" 
	unless ($checkR == 0);
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

## preprocess the simDat
if (!defined($opt_n)) {
    ## find the number of columns
    my $numColInfile = ColNumTabDelimFile($simDat);
    
    # getting column structures from R
    my ($arrRef1, $arrRef2) = GetPriorSumStatNames();
    my @priorNames = @$arrRef1;
    my @sumStatNames = @$arrRef2;
    my $numPriorCols = scalar(@priorNames);

    if (! defined($opt_a)) {  # use the external acceptRejection C program
	my $tol = (defined($opt_t)) ?  $opt_t :  $defaultTolerance;
	if (($numColInfile - @priorNames) % @sumStatNames != 0) {
	  die "ERROR: Simulation file contains $numColInfile columns.\n" .
	      "       R scripts says " . scalar(@priorNames) . " priors and ".
	      scalar(@sumStatNames) . " summary stats for each taxon pairs.\n".
	      "       $numColInfile - " . scalar(@priorNames) . " = " .
	      ($numColInfile - scalar(@priorNames)) . 
	      " should be multiple of ". scalar(@sumStatNames) . ".\n";
	}
	my $numTaxonPairs = ($numColInfile-@priorNames)/scalar(@sumStatNames);

	my @usedSS = split /\s*,\s*/, $statString;
	my @index=();
	
	# finding the column numbers to use as the summary statistics
	my @indexHelper = 1..$numTaxonPairs;
	foreach my $ss (@usedSS) {
	    my $ssi = FindMatchingIndex($ss, @sumStatNames); # 0-offset
	    my @tmp = map { $_ * ($ssi+1)+scalar(@priorNames)} @indexHelper;
	    push @index, @tmp;  # @index is 1-offset
	}
	
	# run rejection program
	my $columns = join " ", @index;
	my $rejExe = FindFile($rejectionExe);
	if ($rejExe eq '-1') {
	    die  "ERROR: Can't find $rejectionExe, is it installed in your PATH?\n";
	}
	system ("$rejExe $obsDat $simDat $tol $columns > $tmpSimDat");

	## create the column only file
	open SIMDAT, "<$simDat" || die "Can't open $simDat\n";
	while (<SIMDAT>) {
	  chomp;
	  my @line = split /\t/;
	  my @priors = splice @line, 0, $numPriorCols;
	  print $tmpPriorfh join("\t", @priors), "\n";
	}
	close SIMDAT;
    }

    my @critVals = (0.01, 0.05, 0.1);
    my @tmpCritVals =  ();
    foreach my $cc (@critVals) {
	push  @tmpCritVals, ($cc) x $numPriorCols;
    }

    # return the proportion of values below the threshold
    # used to calculate bayes factor
    my @priorLT = 
	FreqOfValuesLessThan([1..$numPriorCols],\@tmpCritVals, $simDat);
}

# run R
my @output = `R --quiet --no-save --no-restore --slave < $tmpR`;

print @output;

exit (0);

# find all required R scripts, the main R script has to source several
# R scripts within.  This function will return a text string of the R main script after
# a modification so that source() points to the correct file
sub ProperMainRScript {
    my $result = "";
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

    $result = "source(\"$make_pd\")\nsource(\"$loc2plot\")\n";

    open MAIN_R_TMPL, "<$mainR";
    while(<MAIN_R_TMPL>) {
	s/^\s*source\s*\(\s*["']$make_pdRscript['"]\s*\)\s*$//;
	s/^\s*source\s*\(\s*["']$loc2plotRscript['"]\s*\)\s*$//;
	$result .= $_;
    }
    close MAIN_R_TMPL;

    return $result;
}

# Making R script to print out stat names
sub MkPrintNameRScript {
    my $fh = shift;

    my $mainRScript = ProperMainRScript();
    print $fh "$mainRScript\n";
    if (defined($opt_c)) {
	print $fh "printStatNames(constrained=T)\n";
    } else {
	print $fh "printStatNames(constrained=F)\n";
    }
    return;
}

# making the R script
sub MkStdAnalysisRScript {
    my $fh = shift;

    my $mainRScript = ProperMainRScript();
    print $fh "$mainRScript\n";
    
    if(defined($opt_a)) {
      print $fh "res <- stdAnalysis(\"$tmpObs\", \"$simDat\", pdf.outfile=\"$pdfOut\",pre.rejected=F";

    } else {
      print $fh "res <- stdAnalysis(\"$tmpObs\", \"$tmpSimDat\", \"$tmpPrior\",pdf.outfile=\"$pdfOut\",pre.rejected=T";
    }
    if (defined($opt_a)) {
      my $tol = (defined($opt_t)) ? $opt_t : $defaultTolerance;
      print $fh ", tol=$tol";
    } else {
      print $fh ", tol=1";
    }
    
    print $fh ", used.stats=c($statString)";

    if (defined($opt_r)) {
	print $fh ", rejmethod=T";  # no regression
    } else {
	print $fh ", rejmethod=F";  # regression
    }

    if (defined($opt_c)) {
	print $fh ", constrained=T";  # constrained
    } else {
	print $fh ", constrained=F";  # not constrained
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
	    print STDERR "FILEINFO: using $dir/$name\n";
	    return "$dir/$name";
	}
    }

    return -1;
}

## Read in a tab-delimited text file.
## Then it returns the array whose elements are the frequencies
## of values in columns less than the critical values.
## The column numbers are specified by an array, and the reference to
## the array should be passed as the first argument.
## 2nd argument is the reference to the array of critical values.
## If the two arrays of column numbers and critical values have different
## length, the shorter ones are recycled in the way similar to R.
## Example:
##   FreqOfValuesLessThan([1,3],[0.5,0.5,1,1], "infile.txt")
## returns an array of
##   (freqs of values less than 0.5 in 1st column,
##    freqs of values less than 0.5 in 2nd column,
##    freqs of values less than 1   in 1st column,
##    freqs of values less than 1   in 2nd column)

sub FreqOfValuesLessThan {
    my($colNumArrRef, $valueArrRef, $filename) = @_;

    open IN, "<$filename" || die "Can't open $filename\n";

    my $numCol = scalar(@$colNumArrRef);
    my $numCritVal = scalar(@$valueArrRef);

    # implement R style index recycling
    my @colNums = ();
    my @critVals = ();
    if ($numCol < $numCritVal) {
        for (my $i = 0; $i < $numCritVal; $i++) {
            push @colNums, $$colNumArrRef[$i % $numCol];
        }
        @critVals = @$valueArrRef;
        $numCol = $numCritVal;
    } elsif ($numCritVal < $numCol) {
        $numCol = $numCol;
        @colNums = @$colNumArrRef;
        for (my $i = 0; $i < $numCol; $i++) {
            push @critVals, $$valueArrRef[$i % $numCritVal];
        }
    } else {
        @colNums = @$colNumArrRef;
        @critVals = @$valueArrRef;
    }

    my @result = map {0} 1..$numCol;

    my $cntr = 0;
    while(<IN>) {
        chomp;
        my @line = split /\t/;
        for( my $i = 0; $i < $numCol; $i++) {
            my $index = $colNums[$i] - 1;  # change to 0-offset
            if ($line[$index] < $critVals[$i]) {
                $result[$i] ++;
            }
        }
        $cntr ++;
    }
    close (IN);
    return map {$_ / $cntr} @result;
}

# find the names of Prior columns and sumstats
sub GetPriorSumStatNames {
  my @priorNames = ();
  my @sumStatNames = ();

  my $additionalOpt =  (defined($opt_c))? "-c" : "";
  open NAMES, "$0 -n $additionalOpt 2> /dev/null |" || die "Can't run $0 -n\n";
                                        # running itself to get stat.names.
  my $state = 1;
  while(<NAMES>) {
    chomp;
    next if /^\s+$/;
    if ($state == 1 && /params.from.priorDistn/) {
      $state++; next;
    }
    if ($state == 2) {
      @priorNames = split /\s+/;
      $state++; next;
    }
    if ($state == 3 && /summary.stat.names/) {
      $state++; next;
    }
    if($state == 4) {
      @sumStatNames = split /\s+/;
      last;
    }
  }
  close NAMES;
  return (\@priorNames, \@sumStatNames);
}

# Find the first incidence of target value from the array and return the index
#  Argument: ($target, @array)
sub FindMatchingIndex {
    my $target = shift;
    for (my $i = 0; $i < @_; $i++) {
	return $i if ($target eq $_[$i]);
    }
}

# Find the number of columns in tab delimited text file.
sub ColNumTabDelimFile {
  my $simDat = shift;
  my ($line, $numColInFile);
  open IN, "<$simDat" || die "Can't open $simDat\n";
  while (defined ($line = <IN>)) {
    chomp $line;
    next if $line =~ /^\s*$/;
    my @a = split /\t/, $line;
    $numColInfile = @a;
    last;
  }
  close IN;
  return $numColInfile;
}
