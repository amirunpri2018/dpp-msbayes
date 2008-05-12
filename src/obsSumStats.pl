#! /usr/bin/perl -w

# obsSumStats.pl
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

# version: 20060615
# Originally written by ????, and completely rewritten by 
# Mike Hickerson <mhick@berkeley.edu> and Naoki Takebayashi <ffnt@uaf.edu>
# Licensing: GPL

## dummy strings (DUM_...) will be replaced.
# simpler version.  The most numbers used in this equation doesn't
# influence the output.  Thre of them (-t, 12.0 and 3.1) influence the
# 3rd column of summary stat vector, but this value is not used in later 
# analysis.

my $headerTmpl = <<HEADER_TEMPLATE;
./msDQH 1 DUM_TotSampleSize 1 -t 4.0 -Q 2.0 0.25 0.25 0.25 0.25 -H 999.0 -r 0 DUM_SeqLen -D 5 2 DUM_SampleSize1 DUM_SampleSize2 0 I 0.0 1.0 1.0 0.9 0.9 3.1 2 1 0 0 1 0 I 0.0 Nc 0.4 0.09 12.0 1 Nc 0.03 1 1 Nc 0.03 DUM_SeqLen 1 Nc 0.03 DUM_TaxonPairID 1 Nc 0.03 DUM_NumTaxonPairs

777
11
//
segsites: DUM_NumSeqSites
freqACGT: 0.25 0.25 0.25 0.25
positions:DUM_positions
HEADER_TEMPLATE
###### END of HERE doc ######

my $usage = "Usage: $0 [-h] [-t headerTmplFile] SampleSize_MuModel_Vector\n".
    "  -h: help, print this message\n" .
#    "  -g: Do not remove the sites with any gaps. This probably cause problems.\n".
#    "      So do not use this option." .
    "  -t: an example output file of msDQH (which is used in msbayes.pl)\n" .
    "      You probably do not need to use this option\n";

# Note that the aligned sequence data are represented by a simple 1 dimensional
# array.  Each element is a tab-separated (or a character specified by $sep)
# string, which is in the form of "seqName\tseqData".
# It's more efficient to use array of hash or array of array to represent them,
# but it works.

use strict;
use IO::File;
use POSIX qw(tmpnam);
use Getopt::Std;

# Adding the following paths to @INC, so we can find the R scripts.
# The R scripts should be in the same directory as this perl script,
# or inside of ../lib/msbayes/ relative to this perl script.
# i.e. If this script is inside of ~/bin/, the 3 r-scripts should be inside of
# ~/lib/msbayes/.
use FindBin qw($Bin);
use lib "$Bin";
use lib "$Bin/../lib/msbayes";
use lib ".";  # looks for required files in the current directory at first

my $sumStatsBinName  = "sumstatsvector";
my $sep = "\t";  # used as the internal separater between seq and seq name
my $filename = "SampleSize_MuModel_Vector";  # default filename

our($opt_h, $opt_t, $opt_g);
my ($verbose);  # assign 1 for verbose (for debug).

getopts('hgt:') || die "$usage\n";
die "$usage\n" if (defined($opt_h));


my $sumStatsBin = FindFile($sumStatsBinName);  # locate the binary
if ($sumStatsBin eq '-1') {
    die "Can't find $sumStatsBinName in directories:\n", join(":", @INC), "\n";
}

if (@ARGV > 0) {
    if (@ARGV == 1) {
	$filename = shift;
    } else {
	die "Please give only 1 argument\n";
    }
}

my @obsSumStats = CreateObsSumStats($filename, $headerTmpl);

print join "\n", @obsSumStats;
    
exit;

### This is the actual main function
sub CreateObsSumStats {
    my ($fileName, $header) = @_;

    my @master_matrix = ReadInMaster($fileName);
    my $numTaxonPairs = @master_matrix;
    # create a matrix from the masterinfile (=$fileName) data.
    # In this master matrix following info will be used later:
    #	v DUM_TotSampleSize = column 0 # total sample
    #	w DUM_SampleSize1 = column 1   # taxon A sample
    #	x DUM_SampleSize2 = column 2   # taxon B sample
    #	y DUM_SeqLen = colomn 4        seqLen
    #	u = column 8                   taxon name
    # column numbers here are 0-offset
    
    # Create an array of filenames in the current working directory
    my @directory_files = glob("*");

    # get the header from a file
    # probably we don't need this
    if (defined($opt_t)) {
	$header = GetFileContentsAsScalar($opt_t);	
    }
    
    ## set a temporary file, used as a input for msSumStats
    my ($fh, $tmpFile);
    do {$tmpFile = tmpnam()} 
    until $fh = IO::File->new($tmpFile, O_RDWR|O_CREAT|O_EXCL);
    $fh->close;
    # when we exit or die, automatically delete this temp file
    END{if (defined($tmpFile) && -e $tmpFile) { unlink($tmpFile) || die "Couldn't unlink $tmpFile : $!"}}
    #print "TMP file = $tmpFile\n";  #debug

    my ($scratchfh, $scratchFile);
    do {$scratchFile = tmpnam()} 
    until $scratchfh = IO::File->new($scratchFile, O_RDWR|O_CREAT|O_EXCL);
    $scratchfh->close();
    # when we exit or die, automatically delete this temp file
    END{if (defined($scratchFile) && -e $scratchFile) { unlink($scratchFile) || die "Couldn't unlink $scratchFile : $!"}}
    
    # Each row in @master_matrix corresponds to an alignment file
    # For each row, modify the header template with appropriate values of
    # @master_matrix, and write to $tmpFile.  This file get fed to
    # sumstats program, and receive the results with @sumStatsResultArr.
    my @sumStatsResultArr = ();
     for (my $i = 0; $i < @master_matrix; ++ $i) {
	 ## the following 5 columns are relevant
	 my ($totSampleSize, $sampleSize1, $sampleSize2) = 
	     @{$master_matrix[$i]}[0..2];
	 my ($seqLen, $taxonName) = @{$master_matrix[$i]}[4,8];
	 
	 ### read in the aligned sequence file, and process it.
	 # column 9 (index 8) contains a name for a taxon-pair
	 my $fileName = FindSeqFile($taxonName, \@directory_files);
	 if ($fileName eq "") {
	     die "ERROR: Couldn't open a file for taxon, $taxonName\n";
	 }
	 
	 my @alignedSeq =  (IsFASTA($fileName)) ? ReadInFASTA($fileName) : 
	     ReadInTxt($fileName);
	 
	 # '?' attached to shorter seqs
	 @alignedSeq = AdjustSeqLength(@alignedSeq);
	 
	 unless  (defined($opt_g)) {
	     @alignedSeq = RemoveSitesWithGaps(\@alignedSeq);
	 }
	 @alignedSeq = ExtractVarSites(\@alignedSeq);
	 @alignedSeq = GetSeqDat(@alignedSeq); # get rid of the sequence names
	 
	 ### prepare the header file for this taxon pair
	 my $serialNum = $i + 1;
	 my $numSeqSites = length($alignedSeq[0]);  # number of variable sites
	 
	 my $positionString = "";
	 foreach my $srNum (1..$numSeqSites) {
	     $positionString .= sprintf("%7d", $srNum);
	 }
	 
	 my $header_interim = $header;
	 $header_interim =~ s/DUM_TotSampleSize/$totSampleSize/g;
	 $header_interim =~ s/DUM_SampleSize1/$sampleSize1/g;
	 $header_interim =~ s/DUM_SampleSize2/$sampleSize2/g;
	 $header_interim =~ s/DUM_SeqLen/$seqLen/g;
	 $header_interim =~ s/DUM_NumSeqSites/$numSeqSites/g;
	 $header_interim =~ s/DUM_TaxonPairID/$serialNum/g;
	 $header_interim =~ s/DUM_NumTaxonPairs/$numTaxonPairs/g;
	 # Mike said the next one doesn't matter, but I'm doing it anyway.
	 # It is supposed to contain site positions of mutations
	 $header_interim =~ s/DUM_positions/$positionString/g;
	 
	 ### some consistency check
	 if (@alignedSeq != $totSampleSize) {
	     die "ERROR: taxon, $taxonName, should have " .
		 "$totSampleSize samples.\nBut ", scalar(@alignedSeq),
		 " samples are in the file $fileName\n";
	 }
	 if ($totSampleSize != $sampleSize1+$sampleSize2){
	     die "ERROR: Total sample size ($totSampleSize) for taxon $taxonName".
		 "should be \n       the sum of sampleSizes for the pairs: ".
		 "$sampleSize1 + $sampleSize2.\n".
		 "Check the sampleSize/mutModel File\n";
	 }
	 if ($seqLen < $numSeqSites) {
	     die "ERROR: For taxon, $taxonName, lengths of sequences ($seqLen) " .
		 "should be\n       longer than number of variable sites " .
		 "($numSeqSites)\n";
	 }

	 ### write temporary file as input for sumstat
	 # Basically this file contains the header (fake msDQH command
	 # line) with appropriate values.  After the header, extracted
	 # variable sites data are attached.  It immitates the output
	 # of msDQH
	 open (TMP, ">$tmpFile") || die "Can't open the temp file\n";
	 print TMP $header_interim, "\n", join("\n", @alignedSeq), "\n";
	 close (TMP);

	 warn "INFO: taxon= $taxonName\tfile= $fileName\t# variable sites" .
	     "= $numSeqSites\n";
	 
	 ### run sumstat.
	 my $sumStatsResults = `$sumStatsBin -H --tempFile $scratchFile < $tmpFile`;
	 
	 if ($sumStatsResults !~ /^\s*$/) {  # ignoring empty returned results
	     push @sumStatsResultArr, $sumStatsResults;
	 }
     }
    
    warn "INFO: Number of taxon pairs in the data set = $numTaxonPairs pairs\n";
    return @sumStatsResultArr;
}

# subroutine to get header template in a file and put it into a scalar
# Takes 1 argument: the file name
sub GetFileContentsAsScalar {    
    my $filename = shift;
    my $contents = "";  # scalar with data
    
    open(INPUT, $filename) || die "Cannot open file \"$filename\"!\n";
    while(<INPUT>) {
	$contents .= $_;
    }
    close INPUT;
    return $contents;
}

# Takes a filename as an argument
# Parse the file, and returun a 2-dim matrix
# The file contains information about sample sizes and mutational models.
# It should be tab delimited text file with following columns.
#   1: TotalSampleSize
#   2: SampleSize1 
#   3: SampleSize2
#   4: transition/transversion Ratio
#   5: gamma parameter (not used)
#   6: baseTotalpairs (length of sequences)
#   7: Afreq 
#   8: Cfreq
#   9: Gfreq
#  10: TaxonPairName
# Each line contains a data for 1 taxon pair (sp.1 and sp.2).
# The returned 2-dim matrix contain these information, each line = each row.
# '#' is used to indicate comments, and ignored.
# Empty lines are ignored.
# The configuration info for msprior can be included before this data matrix.
# They include '=', so the first non-empty line which does not contain the
# '=' is considered as the beginning of the samplesize/mut. model data.
sub ReadInMaster {
    my $filename = shift;
    open(INPUT, $filename) || die "Cannot open file \"$filename\"!\n";
    my @master_matrix;
    my $numCol = -1;
    my $paraConfSect = 1;
    my $warnFlag = 1;
    while (<INPUT>) {
	s/#.*$//;  # remove any thing after "#", which indicates comments
	next if (/^\s*$/);
	if ($warnFlag) {   # Maybe we should handle line ending better at some point
	    if (/\r\n$/) {
		warn "WARN: the file $filename appear to have DOS style line ending, it may not work correctly if it's not converted to unix style newline characters\n";
	    }
	    $warnFlag = 0;
	}
	# There could be parameter config lines for msprior in the
	# beginning of this file, which should be ignored.  These are
	# in the format of "parameterName = value", so they are
	# indicated by existence of '=' character. The first line
	# which doesn't contain '=' is considered to be the 1st line
	# with mutation parameter data.
	if ($paraConfSect && /=/) {
	    next;
	} else { # end of msprior para conf section, so lower the flag
	    $paraConfSect = 0;
	}
	
	# when reached here, we are getting into sample sizes/mutational para
	my @master = split /\t/; 
	chomp ($master[$#master]); # remove newline in the taxon name cell
	# check all rows have good column numbers
	if ($numCol < 0) {
	    $numCol = @master;
	    if ($numCol >9 || $numCol < 8) {
		die "ERROR: reading $filename, the 1st line of sample " .
		    "sizes/mutation parameter lines should have 8 or 9 " .
		    "columns.  But, it has $numCol:  $_\n";
	    }
	} elsif ($numCol != @master) {
	    die "ERROR: the number of columns in $filename should be $numCol, " .
		"but the following line has ", scalar(@master), " columns: $_\n";
	}
	
	# total sample numbers should be calculated from the sample sizes of
        # taxon pairs.
	if ($numCol == 8) {
	    unshift @master, $master[0] + $master[1];
	}
	push @master_matrix, [ @master ];
    }
    close INPUT;
    
    return @master_matrix;
}

# This have a problem (original method)
# E.g., if "taxa1.txt" and "taxa12.txt" exist "taxa1" matches with both file
sub FindSeqFileVague {
    my ($name, $arrRef) = @_;
    
    my @matching = grep { $_ =~ /$name/ } @$arrRef;

    if (@matching == 1) {
	return $matching[0];
    } else {
	return "";
    }
}

# Takes the name for the taxon pair as the first arg. and array ref to
# an array containing the name of files in the directory as the 2nd arg.
# It returns the candidate file for the taxon pair.
# If it can't identify a single file, it returns "".
#
# I'm assuming that the filename is taxa1.txt, taxa12.fasta, taxa3,
# etc.  It takes the part of filename before the 1st '.' and check if
# it matches with the name for the taxon pair.  Also exact match is
# ok.  case-insensitive match
sub FindSeqFile {
    my ($name, $arrRef) = @_;
    
    my @matching = grep { ($_ =~ /^$name\./ || $_ =~ /^$name$/) 
			      && $_ !~ /~$/ } @$arrRef;

    if (@matching == 1) {
	return $matching[0];
    } elsif (@matching > 1) {
	warn "ERROR: trying to find the sequence file for taxon pair, $name.\n".
	    "       Multiple matches: ", join(" ", @matching), "\n";
	return "";
    }

    # no match, so ignore the case and see if some file matches
    @matching = grep { $_ =~ /^$name\./i || $_ =~ /^$name$/i} @$arrRef;
    if (@matching == 1) {
	if ($verbose) {
	    warn "Trying to find filename containing $name, similar to " .
		"$matching[0].  This is used\n";
	}
	return $matching[0];
    } elsif (@matching > 1) {
	warn "ERROR: trying to find the sequence file for taxon pair, $name.\n".
	    "       Multiple matches: ", join(" ", @matching), "\n";
	return "";
    }
    # we could use the original relaxed match here, but I'm not doing it
    return "";
}

# each line contains one sequence (sample), no names are give.
sub ReadInTxt {
    my $infile = shift;
    open(INFILE, "<$infile")  || die "Can't open $infile\n";
    
    my @result = ();
    my $baseName = "seq";
    my $cntr = 1;
    while(<INFILE>) {
	chomp;
	next if (/^\s*$/);
	push @result, "$baseName$cntr\t$_";
	$cntr++;
    }
    return @result;
}

# Take a filename as the arg.
# Check if it appear to be fasta file or not.
# Returns:
#   1: the file is FASTA
#   0: the file is not FASTA
# criterion:
#  First non empty line starts with \s*>someName  &&
#  The non empty line after '>' does not conatin '>'
sub IsFASTA {
    my $name = shift;
    open(INFILE, "<$name") || die "Can't open $name\n";
    my $state = 0;
    while(<INFILE>) {
	chomp;
	next if (/^\s*$/);
	
	if ($state == 0) {
	    if (/^\s*>.+$/) {
		$state = 1;
	    } else {
		close INFILE;
		return 0;
	    }
	} elsif ($state == 1) {
	    if (/^\s*>/) {
		close INFILE;
		return 0;
	    } else {
		close INFILE;
		return 1;
	    }
	}
    }
    close INFILE;
    return 0;
}

# takes an arg; name of a file from which data are read Then read in
# the data and make an array.  Each element of this array corresponds
# to a sequence, name tab data.
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";
    
    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
            s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            tr/uU/tT/;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . $_;
        }

        # checking no occurence of internal separator $sep.
        die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
             "the input FASTA file contains this charcter. Make sure this " .
             "separator character is not used in your data file or modify " .
             "variable \$sep in this script to some other character.\n")
            if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
        $result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}

# take std seq data array (name\tseq), and attach "?" at the end for
# shorter sequences
sub AdjustSeqLength {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);

    foreach my $i (0 .. $#seqDat) {
        my $thisLen = length($seqDat[$i]);
        if ($thisLen == $maxLen)  {
            ; # do nothing
        } elsif ($thisLen < $maxLen) {
            my $diff = $maxLen - $thisLen;
            warn "WARN: $seqName[$i] shorter.  " .
                "$diff '?' (missing character) were added at the end\n";
            for (my $j=0; $j < $diff; $j++) {
                $data[$i] = $data[$i] . "?";
            }
        } else {
            die "ERROR: the length of sequence $seqName[$i] is $thisLen, " .
                "longer than \$maxLen = $maxLen.  Weird!!";
        }
    }
    return (@data);
}

# '-' or '?' are considered as gaps
sub RemoveSitesWithGaps {
    my $datArrRef = shift;
    my @seqDat = GetSeqDat(@$datArrRef);
    my @seqName = GetSeqName(@$datArrRef);
    my $maxLen = MaxSeqLen(@$datArrRef);
    my @gapSites = ();
    my @notGapSites = ();
    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
        my @tmpArray = split(//, $seqDat[$seqNumber]);
        # Check the length
        if (@tmpArray != $maxLen)  {
            die "ERROR: the sequence $seqName[$seqNumber] is not same length ".
                "as \$maxLen = $maxLen.  Weird!!";
        }
        push @seqMat, [ @tmpArray ];
    }
    
    # now identify the sites with any gap.
    for $posi (0 .. ($maxLen-1)) {
        my $gap = 0;
        for $seqNumber (0 .. $#seqMat){
            if ($seqMat[$seqNumber][$posi] =~ /^[-\?]$/) {
                $gap = 1;
                last;
            }
        }
        if ($gap == 1) {  #  a gap at these sites
            push (@gapSites, $posi+1); # now unit-offset
        } else {          # there are some non-gap character at these sites
            push (@notGapSites, $posi);
        }
    }
    
    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
        my @thisSeq = SelectSites($seqMat[$seqNumber], \@notGapSites);
        my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
        push (@result, $line);
    }
    
    if ($verbose && @gapSites > 0) {
        warn ("Following sites contains gap(s), removed from analysis\n");
        print STDERR join(" ", @gapSites);
        print STDERR "\n";
    }
    return (@result);
}

## select sites with any variation.  This include varaition due to gaps or 
# degenerate (ambiguous) sites
sub ExtractVarSites {
    my $datArrRef = shift;
    my @seqDat = GetSeqDat(@$datArrRef);
    my @seqName = GetSeqName(@$datArrRef);
    my $maxLen = MaxSeqLen(@$datArrRef);

    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
        my @tmpArray = split(//, $seqDat[$seqNumber]);
        # Check the length
        if (@tmpArray != $maxLen)  {
            die "ERROR: the sequence $seqName[$seqNumber] is not same length ".
                "as \$maxLen = $maxLen.  Weird!!";
        }
        push @seqMat, [ @tmpArray ];
    }
    
    # now identify the sites with variable sites
    my @varSites = ();
    for $posi (0 .. ($maxLen-1)) {
        my $variable = 0;
	my $char1stSeq = $seqMat[0][$posi];
        for $seqNumber (1 .. $#seqMat){
            if ($seqMat[$seqNumber][$posi] !~ /$char1stSeq/i) {
                $variable = 1;
                last;
            }
        }
        if ($variable == 1) {  #  variable site
            push (@varSites, $posi);
        }
    }
    
    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
        my @thisSeq = SelectSites($seqMat[$seqNumber], \@varSites);
        my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
        push (@result, $line);
    }
    
    return (@result);
}

sub SelectSites {
    my ($arrayRef, $indexRef) = @_;
    unless (@_ == 2 && ref($arrayRef) eq 'ARRAY' && ref($indexRef) eq 'ARRAY'){
        die "args to SelectSites() should be ARRAY REF, ARRAY REF\n";
    }

    my $maxIndex = @$arrayRef -1;
    my @result = ();
    foreach my $posi (@$indexRef) {
        if ($maxIndex < $posi) {
            push @result, "?";
        } else {
            push @result, $$arrayRef[$posi];
        }
    }
    return @result;
}

sub MaxSeqLen {
    my @data = GetSeqDat(@_);
    my $maxLen = 0;
    foreach my $i (@data) {
        my $len = length($i);
        $maxLen = $len if ($len > $maxLen);
    }
    return ($maxLen);
}

sub GetSeqDat {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
        @line = split (/$sep/, $i);
        push @result, $line[1];
    }

    return (@result)
}

sub GetSeqName {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
        @line = split (/$sep/, $i);
        push @result, $line[0];
    }
    return (@result)
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
