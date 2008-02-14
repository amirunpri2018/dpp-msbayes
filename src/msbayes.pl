#!/usr/bin/perl

# msbayes.pl
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

my $usage="Usage: $0 [-hd] [-s seed]\n".
    "  -h: help\n".
    "  -d: debug (msprior and msDQH uses the same initial seed = 1)\n".
    "  -s: set the initial seed (but not verbose like -d)\n" .
    "      By default (without -s), unique seed is automaically set from time\n".
    "  -o: output file name\n" .
    "  -c: configuration file for msprior.  Parameters setup interactively,\n".
    "      if this option is not specified\n" .
    "  -r: number of repetitions";


my $defaultOutFile = "Prior_SumStat_Outfile";

use File::Copy;
use IO::File;
use POSIX qw(tmpnam);

use Getopt::Std;

getopts('hdo:c:r:s:') || die "$udage\n";
die "$usage\n" if (defined($opt_h));

my $debug=0;
if (defined($opt_d)) {
    $debug = 1;
}

my $options = "";
if($debug) {  # force msprior to use the same seed
    $options = "-d 1 ";
}

if (defined($opt_r)) {
    if ($opt_r < 0) {
	die "The number of repetitions should be positive integer\n$udage\n";
    }
    $options = $options . " --reps $opt_r ";
}

if (defined($opt_c)) {
    die "ERROR: $opt_c is not readable\n" unless (-r $opt_c);
    die "ERROR: $opt_c is empty\n" if (-z $opt_c);
    die "ERROR: $opt_c is not a text file\n" unless (-T $opt_c);
    $options = $options . " --config $opt_c ";
}

if (defined ($opt_s)) {
    $options = $options . " --seed $opt_s ";
}

if (defined($opt_s) || defined($opt_d)) {  # set the msDQH use the same seeds
    if (defined($opt_s)) {
	srand($opt_s);
    } else {
	srand(1);
    }
}

#### Find programs
my $msprior = FindExec("msprior");
my $msDQH = FindExec("msDQH");
my $sumstatsvector = FindExec("sumstatsvector");

# open and close a temp file
# This is used to store the prior paras from msprior (psiarray and tauarray)
my $tmpPriorOut, $tmpPriorOutfh;
do {$tmpPriorOut = tmpnam()} until $tmpPriorOutfh = 
    IO::File->new($tmpPriorOut, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpPriorOut) && -e $tmpPriorOut) {
	unlink($tmpPriorOut) || die "Couldn't unlink $tmpPriorOut : $!"
	}
};
$tmpPriorOutfh->close();
$options = $options . " --priorOut $tmpPriorOut ";

# The main results (after sumstatsvector) get stored in this file
# Then this and $tmpPriorOut files get column concatenated produce the final
# output.  As long as the /tmp is local file, this enable running the
# program in NFS mounted /home
my ($tmpMainOut, $tmpMainOutfh);
do {$tmpMainOut = tmpnam()} until $tmpMainOutfh = 
    IO::File->new($tmpMainOut, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpMainOut) && -e $tmpMainOut) {
	unlink($tmpMainOut) || die "Couldn't unlink $tmpMainOut : $!"
	}
};
$tmpMainOutfh->close();

# This is used by sumstats to store temporary output.  It used to be
# "PARarray-E", but now we are using temp file.  Better for NFS home
my ($tmpSumStatVectScratch, $tmpSumStatVectScratchFh);
do {$tmpSumStatVectScratch = tmpnam()} until $tmpSumStatVectScratchFh = 
    IO::File->new($tmpSumStatVectScratch, O_RDWR|O_CREAT|O_EXCL);
END {                   # delete the temp file when done
    if (defined($tmpSumStatVectScratch) && -e $tmpSumStatVectScratch) {
	unlink($tmpSumStatVectScratch) || die "Couldn't unlink $tmpSumStatVectScratch : $!"
	}
};
$tmpSumStatVectScratchFh->close();

#### setting up output filename
my $outFile;
if(defined($opt_o)) {
    $outFile = $opt_o;
    CheckNBackupFile($outFile);
} else {
    $outFile = InteractiveSetup();
}

my $new = 1;
open (RAND, "$msprior $options |") || 
    die "Hey buddy, where the hell is \"msprior\"?\n";

## Note this file is used as temp file in sumstats.  We need to clean it
## before computation if there is a leftover.
CheckNBackupFile("PARarray-E");

my $counter  = 1;
while (<RAND>) {
    s/^\s+//; s/\s+$//; 

    my ($upperTheta, $theta, $gaussTime, $mig, $rec, $taxonPairID,
	$BottleTime, $BottStr1, $BottStr2, 
	$totSampleNum, $sampleNum1, $sampleNum2, $tstv1, $tstv2, $gamma,
	$seqLen, $numTauClasses, $N1, $N2, $Nanc, 
	$freqA, $freqC, $freqG, $freqT, $numTaxaPair) = split /\s+/;

# 0 $upperTheta, $theta, $gaussTime, $mig, $rec, $taxonPairID,
# 6 $BottleTime, $BottStr1, $BottStr2, $totSampleNum, $sampleNum1,
#11  $sampleNum2, $tstv1, $tstv2, $gamma, $seqLen,
#16  $numTauClasses, $N1, $N2, $Nanc, $freqA,
#21  $freqC, $freqG, $freqT

    my $tmpVal = $gaussTime - $BottleTime;

    # option for -r was fixed to 0, so changed to $rec, then forcing
    # it to be 0 here
    $rec = 0;

    $SEED = int(rand(2**32));  # msDQH expect unsigned long, the max val (2**32-1) is chosen here

    # Printing the header at the right time
    my $headerOpt = ($counter == $numTaxaPair) ? "-H":"";

#    print("$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $tmpVal 1 Nc $Nanc $numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonPairID 1 Nc $Nanc $numTaxaPair | $sumstatsvector -T $upperTheta $headerOpt >> $tmpMainOut"); print ("\n");

    system("$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $tmpVal 1 Nc $Nanc $numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonPairID 1 Nc $Nanc $numTaxaPair | $sumstatsvector -T $upperTheta --tempFile $tmpSumStatVectScratch $headerOpt >> $tmpMainOut");

#   system("$msDQH $SEED $totSampleNum 1 -t $theta -Q $tstv1 $freqA $freqC $freqG $freqT -H $gamma -r $rec $seqLen -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $BottStr1 $N2 $BottStr2 $BottleTime 2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $tmpVal 1 Nc $Nanc $numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonPairID 1 Nc $Nanc $numTaxaPair >> $tmpMainOut");


# minor change 9/8/06; $N1 $N1 $N2 $N2 to $N1 $BottStr1 $N2 $BottStr2

# The command line format is the same as Dick Hudson's ms except for
# everything after -D. Everything after -D specifies what happens
# during X number of time intervals. -D 6 2 means 6 intervals starting
# with 2 populations at time 0 (going backwards in time).

# The rest -D is explained using the following template:
#
# -D 6 2 $sampleNum1 $sampleNum2 0 I $mig $N1 $N1 $N2 $N2 $BottleTime \
#  -2 1 0 0 1 0 I $mig Nc $BottStr1 $BottStr2 $tmpVal 1 Nc $Nanc \
#   -$numTauClasses 1 Nc $Nanc $seqLen 1 Nc $Nanc $taxonPairID 1 Nc
#    -$Nanc $numTaxaPair
#
# $sampleNum1 $sampleNum2; the 2 sample sizes of the 2 populations

# 0 I $mig;  $mig is migration rate under an island model of migration

# $N1 $N1 $N2 $N2; Relative size pop1 at begening of timestep,
#   Relative size pop1 at end of timestep, Relative size pop2 at
#   begening of timestep, Relative size pop2 at end of timestep

# $BottleTime; length of first time step (begining of bottleneck going
#    backwards in time)

# 2 1 0 0 1; this is the admixture matrix a[][]; this allows
#   population divergence or admixture, or can specify neither occuring
#   during the #time step. In this case nothing happens (2 populations
#   become two #populations)
#   1 0
#   0 1
#   in a[][] this means all of pop1 goes into pop1, and all of pop2 goes
#   into pop 2 (2 populations remain isolated)
#   If we had:
#   2 1 0 1 0 0 1
#   1 0
#   1 0
#   0 1
#   this would mean that at first we have 3 populaations and then all of
#   pop2 fuses with pop1 (going back in time and hence divergence), and
#   pop3 would remain intact

# 0 I $mig;  the migration region of the next time step

# Nc; Nc specifies that all populations have constant size in this
#   next time step

# $BottStr1 $BottStr2; these are the two constant realtive sizes of
#   the two populations during this next time step.

# $tmpVal; this is the length of this next time step (in this case it
#   ends at the divergence time)

# 1 Nc $Nanc $numTauClasses; specifies that the next time step has
#   only one population (divergence) and the population is constant in
#   size through the time step
#     $Nanc; relative size of this ancestral population
#     $numTauClasses; this is the length of the time step, but has a 2nd
#        meaning unrelated to msDQH. the actual value gets passed on to the
#        summary stats program for parameter estimation purposes. The actual
#        value is somewhat arbitray for msDQH because there is only one
#        population remaining going back in time.  The length of the period
#        can be infinite.

# Three more time-steps use the same "1 Nc $Nanc $length" pattern,
# where "length" has a 2nd use

# If one wants to use msDQH independently on the command line, one can
# add "-P" to see what population model is being used.  Example below
#
# ./msDQH 35 1 -t 20.0 -Q 5.25 0.25 0.25 0.25 0.25 -H 999.000000 -r 0 1000 -D 5 2 20 15 0 I 0.000000 0.8 0.05 0.9 0.05 6.03 2 1 0 0 1 0 I 0.000000 Nc 0.05 0.05 0.001 1 Nc 0.42 6 1 Nc 0.42 1000 1 Nc 0.42 1 -P
 
# Output example using "-P"

# In the below example 2 populations (20 and 15 individuals) diverged
# from a common ancestor of relative size 0.42 at the third time step

#./msDQH 35 1 -t 20.0 -Q 5.25 0.25 0.25 0.25 0.25 -H 999.000000 -r 0 1000 -D 5 2 20 15 0 I 0.000000 0.8 0.05 0.9 0.05 6.03 2 1 0 0 1 0 I 0.000000 Nc 0.05 0.05 0.001 1 Nc 0.42 6 1 Nc 0.42 1000 1 Nc 0.42 1 -P 
#./msDQH nsam 35 howmany 1
#  theta 20.00 segsites 0
#seQmut 1 output 0
#tstvAG CT 5.25 5.25, freqACGT 0.25 0.25 0.25 0.25
#gammaHet alpha 999.00
#  r 0.00 f 0.00 tr_len 0.00 nsites 1000
#  Dintn 5 
#    Dint0 npops 2
#      config[] 20 15 
#      Mpattern 0
#      M[][] 0.00 0.00 0.00 0.00 
#      (Nrec_Npast)[] 0.80 0.05 0.90 0.05 
#       tpast 6.03
#    Dint1 npops 2
#      a[][] 1.00 0.00 0.00 1.00 
#      Mpattern 0
#      M[][] 0.00 0.00 0.00 0.00 
#      (Nrec_Npast)[] 0.05 0.05 0.05 0.05 
#       tpast 6.03
#    Dint2 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 12.03
#    Dint3 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 1012.03
#    Dint4 npops 1
#      (Nrec_Npast)[] 0.42 0.42 
#       tpast 1013.03
    $counter++;
}

close RAND;

# combine the two outputfile to create the final output file
my $rc = ColCatFiles($tmpPriorOut, $tmpMainOut, $outFile);
if ($rc != 1) {
    if ($rc == 2) {
	warn "WARN: the number of lines in the prior output different from the main".
	    " output.  This should not have happened.  Something went wrong, so ".
	    "DO NOT TRUST THE RESULTS!!!"
    }
    if ($debug && $rc == -1) {    # debug, remove this later
	warn "In Concatenating at the end, one file was empty.  This means that the simulation was the UNCONSTRAINED\n";
    }
}

if (-e "PARarray-E") {
    unlink("PARarray-E") || die "Couldn't unlink PARarray-E : $!";
}
exit(0);

# interactively setting up.
sub InteractiveSetup {
    my $outFileName;

    # output filname
    print "Output Filename? [Return] to use default of " .
	"\"$defaultOutFile\"\n";
    chomp($outFileName = <STDIN>);
    if($outFileName eq "") {
	$outFileName = $defaultOutFile;
    }
    CheckNBackupFile($outFileName);

    return ($outFileName)
}

# This fucntion check if the argument (fileName) exists. If it exists,
# it get renamed to fileName.oldN, where N is a digit.
# In this way, no files will be overwritten.
sub CheckNBackupFile {
    my $fileName = shift;
    if ($fileName eq "") {
	$fileName="Prior_SumStat_Outfile";
    }

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

### Supply the name of program, and it will try to find the executable.
### In addition to regular path, it search for several other places
sub FindExec {
    my $prog = shift;
    # I'm making it to find the binary in the current directory (.) at first.
    # I do not personally like this.  But since Mike doesn't like to
    # install the binaries in the appropriate directories, we need to
    # force this behavior to reduce confusion. 
    # When this program become more matured, we should reevaluate this.
    # Similar behavior in acceptRej.pl introduced  Naoki Feb 8, 2008
    $ENV{'PATH'} = ".:" . $ENV{'PATH'} . 
	":/bin:/usr/bin:/usr/local/bin:$ENV{'HOME'}/bin";
    my $bin = `which $prog 2>/dev/null`;
    chomp $bin;

    if ($bin eq "") {
	die "ERROR: $prog not found in PATH $ENV{'PATH'}\n";
    }

    print STDERR "INFO: using $bin\n";
    return $bin;
}

# Take filenames of two files, and concatenate them side by side to
# produce the outputfile given as the 3rd argument.  The tab will be
# inserted after each line of the first file.
# If two files do not have the same number of lines, the output file
# have the same length as the files with less lines.  The rest of the
# longer file is ignored.

# Return value
# 1 two files same number of lines, and successfully concatenated
# 2 two files have different number of lines
# -1 Either one file or bot files were empty.  non-empty file is copied as the 
#    output file
sub ColCatFiles {
    my ($infilename1, $infilename2, $outfilename) = @_;

    # check empty files, if empty, copy is enough
    if ( -s $infilename1 && -z $infilename2) { # file 2 empty
	copy($infilename1, $outfilename) ||
	    warn "WARN: copy $infilename1, $outfilename failed";
	return -1;
    } elsif (-z $infilename1 ) { # file 1 empty or both empty
	copy($infilename2, $outfilename) ||
	    warn "WARN: copy $infilename2, $outfilename failed";
	return -1;
    }

    # both infiles are not empty
    my $retval = 1;

    open FILE1, "<$infilename1" || die "Can't open $infilename1\n";
    open FILE2, "<$infilename2" || die "Can't open $infilename2\n";
    open OUT, ">$outfilename" || die "Can't open $outfile\n";

    $numLines1 = `wc -l < $infilename1`;
    die "wc failed: $?" if $?;
    chomp $numLines1;
    $numLines2 = `wc -l < $infilename2`;
    die "wc failed: $?" if $?;
    chomp $numLines2;

    if ($numLines1 != $numLines2) {
	warn "WARN: number of lines differ between $infilename1 and $infilename2\n";
	$retval = 2;
    }
    
    my $maxLines =  ($numLines1 > $numLines2) ? $numLines1 : $numLines2;
    
    for(my $i = 0; $i < $maxLines; $i++) {
	if ($i < $numLines1) {
	    $line = <FILE1>;
	    chomp $line;
	    print OUT $line;
	}
	print OUT "\t";

	if ($i < $numLines2) {
	    $line = <FILE2>;
	    chomp $line;
	    print OUT $line;
	}
	print OUT "\n";
    }
    
    close OUT;
    close FILE1;
    close FILE2;
    return $retval;
}
