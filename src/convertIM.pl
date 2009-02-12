#!/usr/bin/perl

#use strict;  # need to check the exact syntax
use File::Copy;
use IO::File;

$fileSuffix;
$theHeader;
$theLine;
$readingHeader = 1;
$readingData = 0;
$startToReadIM = 1;
$batchFileName;
$pop1;
$pop2;


#print "convertIM started \n";
unless (open LISTFILE, "< $ARGV[0]")
{
   die "Can't open file $ARGV[0] : $! \n";    
}   
my @fileList = <LISTFILE>;
close LISTFILE;

foreach my $currentFile (@fileList)
{ 
   # trim the line
   chomp $currentFile;         # remove newline
   $currentFile =~ s/^\s+//;   # remove leading whitespace
   $currentFile =~ s/\s+$//;   # remove trailing whitespace
   
   # check whether currentFile is OK
   die "ERROR: $currentFile is not readable\n" unless (-r $currentFile);
   die "ERROR: $currentFile is empty\n" if (-z $currentFile);
   die "ERROR: $currentFile is not a text file\n" unless (-T $currentFile);
   
   if($startToReadIM == 1)
   {
        # create batch.masterIn file
        $batchFileName = "batch.masterIn.fromIM";

	if (-e $batchFileName)
	{
	    my $i = 1;
	    while (-e "$batchFileName.old$i")
	    {
		  # checking if the file exists
		  $i++;
	    }
	    move ("$batchFileName", "old$i.$batchFileName") ||
		 die "Can't rename $batchFileName to $batchFileName.old$i";
	}

	unless (open BATCH, "> $batchFileName")
	{
	    die "Can't open file $batchFileName : $! \n";    
	}   

	# print mutation table at the begining of batch file
        printf BATCH "upperTheta = 30.0\nupperTau = 10.0\nnumTauClasses = 0\nupperMig = 0.0\nupperRec = 0.0\nupperAncPopSize = 0.5\nnumLoci = 1\nconstrain = 0\nsubParamConstrain = 111111111\n\n";
        printf BATCH ("BEGIN SAMPLE_TBL\n");
        $startToReadIM = 0;
   } # $startToReadIM == 1    


   if( ($fileSuffix = index($currentFile, ".im")) == (length($currentFile)-3) )
   {
          my $pop1Name,$pop2Name,$lineRead = 0,$numLoci,$locus_name, $number_pop1,$number_pop2, $sample_length, 
             $mutation_model,$inheritance_scalar, $ploidy, $line,$A = 0,$C = 0,$T = 0,$G = 0,$total,$allele_count;
              
	  # open a .im file, grap the data, then close it (to save resource)
	  unless (open IMFILE, "< $currentFile")
	   {
	      die "Can't open file $currentFile : $! \n";    
	   }    

	   my @imData =  <IMFILE>;

	   close IMFILE;

	   foreach $line (@imData)  
	   {
		 # trim the line
		 chomp $line;         # remove newline
		 $line =~ s/^\s+//;   # remove leading whitespace
		 $line =~ s/\s+$//;   # remove trailing whitespace

		 # find out whether the line is a comment or an empty line 
		 my $poundFound = index($line, "#");

		 if( $poundFound == 0 && $lineRead != 0)
		 {
		      next;
		 } # the line is a comment

		 # increment $lineRead for each valid line
		 $lineRead++ ;

		 if($lineRead == 1)
		 {
	            next;
		 } # 1st line contains arbitary text

		 if($lineRead == 2)
		 {
	            ($pop1Name, $pop2Name) = split /\s+/, $line;
		     next;    
		 } # 2nd line contains population names

		 if($lineRead == 3)
		 {
		     $numLoci = $line;      
		     next;       
		 } # 3rd line contains number of loci

		 if($lineRead >= 4)
		 {
		    if($readingHeader == 1)
	            {

			 ($locus_name, $number_pop1, $number_pop2, $sample_length, 
			  $mutation_model, $inheritance_scalar) = split /\s+/, $line;
                          
                          if($inheritance_scalar == 1)
			  {
			      $ploidy = 2;
                          }else
                          {
                              $ploidy = 1;
                          }      

			  $total = $number_pop1 + $number_pop2;
 			  $readingHeader = 0;
			  $readingData = 1;

			  next;

		     } # this line contains a header

		     if($readingData == 1)
		     {
                        if($allele_count < $total)
			{
			    my ($tname, $allele) = split /\s+/, $line;

			    my $i = rindex($tname,"_");
			    $name = substr $tname, 0, $i;

			    if ($pop1 eq undef && $pop2 eq undef)
			    {
			       $pop1 = $name;  
			       push(@pop1_data, $allele);
			    }
			    elsif($name eq $pop1)  
			    { 
				push(@pop1_data, $allele);
			    }elsif ($pop2 eq undef)
			    {
				$pop2 = $name;
				 push(@pop2_data, $allele);
			    }elsif($name eq $pop2)
			    {
				 push(@pop2_data, $allele);
			    }

			    # get precentage of A,T,C,G
			    my @charArray = split //, $allele;
			    foreach my $cc (@charArray)
			    {
		              if($cc eq "A")
				   {
				     $A++;
				   }elsif ($cc eq "T"){
				     $T++;
				   }elsif ($cc eq "C"){
				     $C++;
				   }elsif ($cc eq "G"){
				    $G++;
				   }
				  else{
				    # print STDERR "found something weired $cc in an allele sequence \n";
				   }
			     } # foreach in @charArray

			     $allele_count++ ;

			     if($allele_count < $total)
			     {
				  next;
			     }
			  } # $allele_count < $total

			} # $readingData

		    } # $lineRead >= 4

		    # reinitialize some variables
		    $allele_count = 0;
		    $readingData = 0;
		    $readingHeader = 1;

		    # calculate the proportion of A,T,C and G
		    my $totalSampleLength = $sample_length * $total;
		    if ($totalSampleLength == 0)
		    {
			$A = $C = $T = $G = 0;
		    }
		    else
		    {
		     $A = $A / $totalSampleLength;
		     $C = $C / $totalSampleLength;
		     $T = $T / $totalSampleLength;
		     $G = $G / $totalSampleLength;
		    }

		    # construct the header for batch.masterIn file
		    my $ii = rindex($currentFile, ".");
		    my $theAllele = substr $currentFile, 0, $ii;
		    my $fastaFileName = $theAllele . "." . $locus_name . ".fasta";
		    $theHeader = sprintf ("%s\t%s\t%d\t%d\t%d\t1.00\t%d\t%.3f\t%.3f\t%.3f\t%s\n", $theAllele, $locus_name, $ploidy, $number_pop1, $number_pop2, $sample_length, $A, $C, $G, $fastaFileName );
		    printf BATCH ($theHeader);
       		    $A = $C = $T = $G = 0;

		    if (-e $fastaFileName)
		    {
			my $i = 1;
			while (-e "$fastaFileName.old$i")
			{
			  # checking if the file exists
			    $i++;
			}
			move ("$fastaFileName", "old$i.$fastaFileName") ||
			    die "Can't rename $fastaFileName to $fastaFileName.old$i";
		    }

		    unless (open FASTA, "> $fastaFileName")
		    {
			 die "Can't open file $fastaFileName : $! \n";    
		    }   

		    my $theAllele1 = $theAllele . "_" . $pop1Name . "_" . $locus_name;
		    for (my $i = 0; $i < $number_pop1; $i ++)
		    {
		       $theLine = sprintf (">%s\n%s\n", $theAllele1, shift(@pop1_data));
		       printf FASTA ($theLine);
		    }

		    my $theAllele2 = $theAllele . "_" . $pop2Name . "_" . $locus_name;
		    for (my $i = 0; $i < $number_pop2; $i ++)
		    {
		       $theLine = sprintf (">%s\n%s\n", $theAllele2, shift(@pop2_data));
		       printf FASTA ($theLine);
		    }
          
		    close FASTA;

	   } # while (still a line in the current file)
	   
	 $pop1 = $pop2 = undef; 
    } # read in a IM file
    else
    {
       die "$currentFile is not a valid file \n";
    }
            
} # while (still files in the file list)


printf BATCH ("END SAMPLE_TBL\n\n");

printf BATCH ("BEGIN CONSTRAIN\n");

my @csData = ("1.0	0.9	0.1	0.5	0.0	10.1	1.5	0.1	0.0",
 	      "1.1	0.8	0.2	0.6	0.0	20.1	1.4	0.2	0.0",
	      "1.2	0.7	0.3	0.7	0.0	30.1	1.3	0.3	0.0",
	      "1.0	0.3	0.7	0.8	0.0	40.1	1.2	0.4	0.0",
	      "1.0	0.3	0.8	0.9	0.0	5.1	1.1	0.5	0.0",
	      "1.0	0.3	0.9	0.3	0.0	25.1	1.0	0.5	0.0");

foreach my $csLine (@csData)
{
   # trim the line
   chomp $csLine;         # remove newline
   $csLine =~ s/^\s+//;   # remove leading whitespace
   $csLine =~ s/\s+$//;   # remove trailing whitespace
   printf BATCH "$csLine\n";
}
printf BATCH ("END CONSTRAIN\n");

close BATCH;

print "$batchFileName\n";

exit(0);
