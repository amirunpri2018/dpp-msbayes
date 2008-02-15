/*
 *  msprior.c
 *
 * Copyright (C) 2006  Michael Hickerson and Naoki Takebayashi
 *
 * This file is a part of msprior.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


/*
 * This program read in 3 positive integer values as the input.
 *   msprior reps taxa loci
 * reps is number of replications, taxa is the number of taxon pairs
 * loci is number of loci.
 * Additionally, it reads in a configuration file.
 *
 * Output: Each line is tab delimited, and the columns contain following info.
 * 1. theta, from flat prior
 * 2. tau, time when the pair got separated (with weird scaling)
 *         Several classes of tau are chosen from flat prior, and tau for
 *         a species get chosen from these classes.
 * 3. migration rate, fixed value (=0).
 * 4. recombination rate, from flat prior
 * 5. integer index: 1 ... number of taxa
 * 6. BottleTime
 * 7. BottStr1
 * 8. BottStr2
 *
 */

/*
  Change Log
  * Fri May 12 2006 Naoki Takebayashi <ffnt@uaf.edu>
  - upper and lower bounds of prior distribution for theta is configurable
    uses gParam.upperTheta and gParam.lowerTheta
  - To make this to work, tauequlaizer uses gParam.upperTheta/2.  upperTheta
    needs to be passed to the stat program.  So the first element of the 
    output line of this program is upperTheta now.  msbayes.pl correctly
    receive this and pass this value to the option flag (-T) of the 
    stat program.

    * Tue Mar 14 2006 Naoki Takebayashi <ffnt@uaf.edu>
  - separated the function of reading in the config file, and setup parameters.
    This is now done in a separate file setup.c.  Setting up parameters
    are cleaner now.
  - As the consequence, the number of taxon-pairs are not arbitrary limited
  - got rid of unused variables, and general clean-up
  - BottStr1 and 2 were forced to be 0.01 with probabilty 0.1.  This created
    rather weird prior distributions for the strengths of bottleneck.
    This is removed now.
  - Compared the results of the previous and new implementations.  Other than
    this arbitrary preference to stronger bottleneck, the results are 
    completely identical with the same seed.
*/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_rng.h>/* for base rand num gen's */
#include <gsl/gsl_randist.h>/* for gsl_ran_gamma */
#include "msprior.h"
#include "setup.h"

/* This has to be Global */
const gsl_rng *gBaseRand;/* global rand number generator */

/* 
 * global variables which stores the parameters (settings for the upper
 * limist of prior dist'n etc
 */
runParameters gParam;        /* stores upper limits of prior dist'n etc */
mutParameterArray gMutParam; /* stores mut model & # samples for 
                                each taxon pair */

constrainedParameterArray gConParam; /* store constrained sub-parameters for
                                        each taxon pair */ 

int comp_nums(const double *num1, const double *num2)
{
  if (*num1 <  *num2) return -1;
  else if (*num1 > *num2) return 1;
  else return 0;
} 

int main (int argc, char *argv[])
{
  double N1, N2, Nanc, *tauArray=NULL, theta, tauequalizer, gaussTime,
    mig, rec, BottStr1, BottStr2,BottleTime;
  int  tauClass, *PSIarray=NULL;
  unsigned int numTauClasses=-1, u, locus, zzz, c;
  unsigned long randSeed;
  unsigned long long rep;
  extern const gsl_rng *gBaseRand;
  void display_nums(int *, int);
  int comp_nums(const double *, const double *);
  FILE *fpTauPsiArray;

  /* set up gParam and gMutParam */
  LoadConfiguration(argc, argv);

  /* for initiating the gsl random number generator */
  /* initialize PRNG */
  srand (gParam.prngSeed);/* Better way of seeding here ? */
  randSeed = rand ();
  if (debug_level > 0)
    randSeed = 1;

  gBaseRand = gsl_rng_alloc (gsl_rng_mt19937);/* set the base PRNG to
						 Mersenne Twister */
  gsl_rng_set (gBaseRand, randSeed);/* seed the PRNG */

  /* fixed numTauClasses configuration */
  if (gParam.numTauClasses != 0) {
    if (gParam.numTauClasses > gParam.numTaxaPair) {
      fprintf(stderr, "WARN: numTauClasses (%u) is larger than "
	      "numTaxaPair (%u). Setting numTauClasses to %u",
	      gParam.numTauClasses,gParam.numTaxaPair,gParam.numTaxaPair);
      gParam.numTauClasses = gParam.numTaxaPair;
    }
    numTauClasses = gParam.numTauClasses;
  }

  /* Sizes are set to the number of taxon pairs (Max number of tau)*/
  tauArray = calloc(gParam.numTaxaPair, sizeof(double));
  if (tauArray == NULL) {
    fprintf(stderr, "ERROR: Not enough memory for tauArray\n");
    exit(EXIT_FAILURE);
  }
  
  PSIarray = calloc(gParam.numTaxaPair, sizeof(int));
  if (PSIarray == NULL) {
    fprintf(stderr, "ERROR: Not enough memory for PSIarray\n");
    exit(EXIT_FAILURE);
  }

  /* open the stream to output tauArray and PSIarray */
  if ((fpTauPsiArray=fopen(gParam.priorOutFile, "w")) == NULL) {
    fprintf(stderr,"Cannot open the file.\n");
    exit(1);
  }

  /* print out the column headers to the tauArray and PSIarray  */
  /* PRI.numTauClass PRI.Tau.1  PRI.Tau.2 PRI.Tau.3 ... PRI.Psi.1 PRI.Psi.2 PRI.Psi.3 ... */
  fprintf(fpTauPsiArray, "PRI.numTauClass\t");
  if (gParam.numTauClasses > 0) {  /* constrained psi analysis */
    for (zzz = 0; zzz < numTauClasses; zzz++) {
      fprintf(fpTauPsiArray, "PRI.Tau.%d\t",zzz+1);
    }
    for (zzz = 0; zzz < numTauClasses; zzz++) {
      fprintf(fpTauPsiArray, "PRI.Psi.%d",zzz+1);
      fprintf(fpTauPsiArray, ((zzz != numTauClasses - 1) ? "\t" : "\n"));	  
    }
  }
  
  /* Beginning of the main loop */
  for (rep = 0; rep < gParam.reps; rep++)
    {
      /*
       * Each taxon pair was separated at a time tau in the past.  Of
       * all pairs, some of them may have been separated at the same
       * time.  numTauClasses is the number of classes with different
       * divergence time.
       *
       * If gParam.numTauClasses is not set, we are sampling
       * numTauClasses from a uniform prior dist'n.
       */ 
      if (gParam.numTauClasses == 0) { /* numTauClasses is NOT fixed */
	numTauClasses = 1 + gsl_rng_uniform_int(gBaseRand, gParam.numTaxaPair);
      } 
      
      /* sample tau's from uniform prior dist'n */      
      for (u = 0; u < numTauClasses; u++)
	{
	  tauArray[u] = gsl_ran_flat (gBaseRand, 0.0, gParam.upperTau);
          
          //printf("tauArray[%d] : %lf   ", u, tauArray[u]);
          
	  if (debug_level) {
	    fprintf(stderr, "DEBUG:%u of %u categories:\t%lf\n",
		    u, numTauClasses, tauArray[u]);
	  }
	}

      //qsort(tauArray, (numTauClasses), sizeof(double), comp_nums);

      for (c=0; c < numTauClasses; c++) 
	PSIarray[c] = 0;  /* Reset the PSIarray counters */

      for (u = 0; u < gParam.numTaxaPair; u++)
	{
	  mutParameter taxonPairDat;

	  /* sample sizes, mutational model for u-th taxon-pair */
	  taxonPairDat = gMutParam.data[u];

	  /** bottleneck priors **/
	  /* severity of bottle neck (how small the population become) */
	  BottStr1 = gsl_ran_flat (gBaseRand, 0.01, 1.0);
	  BottStr2 = gsl_ran_flat (gBaseRand, 0.01, 1.0);
	  /* timing of bottle neck */
	  BottleTime = gsl_ran_flat (gBaseRand, 0.000001, 1.0);

	  /* migration rate prior */
	  mig = gsl_ran_flat (gBaseRand, 0.0, gParam.upperMig);

	  /* theta prior */
	  theta = gsl_ran_flat (gBaseRand, gParam.lowerTheta, 
				gParam.upperTheta);

	  /* population sizes immediately after the separation, and 
	     what it grows to after the bottleneck (today) */
	  N1 = gsl_ran_flat (gBaseRand, 0.01, 1.99);

	  N2 = 2.0 - N1;
	    
	  /* ancestral population size prior */
	  if(gParam.upperAncPopSize < 0.01) {
	    fprintf(stderr, "The upper bound (%lf) of ancestral pop. size is "
		    "smaller than the lower bound (0.01)\n",
		    gParam.upperAncPopSize);
	    exit(EXIT_FAILURE);
	  }
	    
	  /*
	        Nmax=((gParam.upperAncPopSize*gParam.upperTheta)*gParam.upperTheta)
		    /theta; 
		      Nanc = gsl_ran_flat (gBaseRand, 0.01, Nmax);
	  */
	    
	  /* The upper limit of ancestral theta is defined by the product
	     of upper Theta (e.g. 40) and upper AncPopSize (e.g. 0.5) */
	  Nanc = gsl_ran_flat(gBaseRand, 0.01, 
			      gParam.upperAncPopSize*gParam.upperTheta);
	  Nanc = Nanc / theta; /* get the ratio of theta_anc / theta_cur 
				  This ratio is required for msDQH */

	  tauequalizer = gParam.upperTheta / 2 / theta;
	   
	  /* pick a tau for every taxon-pair with replacement from the
	          array of X taxon-pairs, where X is a uniform discrete RV
		  from 1 to number of taxon-pairs */
	  

          // 1-31-2007 it will use new way of picking tau
	  if( u < numTauClasses)
	    tauClass = u;
	  else
	    tauClass = gsl_rng_uniform_int(gBaseRand,numTauClasses);
	  
	  gaussTime= tauArray[tauClass];
          //printf("picking index: %d, gaussTime: %lf  ", tauClass, tauArray[tauClass]);

	  PSIarray[tauClass] = PSIarray[tauClass] + 1;  

	    
	  /* use the following if simulating a particular fixed history */
	  /* gaussTime = tauArray[u]; */
	    
	  gaussTime = gaussTime * tauequalizer;
	    
	  /* The following 2 if's are weird */
	  if (gaussTime < 0.0001)
	    gaussTime = 0.0001;

	  BottleTime = BottleTime * 0.95 * gaussTime;

	  if (gaussTime < 0.0001)
	    BottleTime = 0.00005;

	  if(debug_level) 
	    fprintf(stderr, "DEBUG: BottleTime:%lf\tgaussTime:%lf\n",
		    BottleTime, gaussTime);
	    
	  /* recombination rate */
	  rec = gsl_ran_flat (gBaseRand, 0.0, gParam.upperRec);
	    
	  /* print out the results */
	  for (locus = 0; locus < gParam.numLoci; locus++) {
	    printf("%lf %lf %lf %lf %lf %u ",
		   gParam.upperTheta, theta, gaussTime, mig, rec, u+1);
	    printf("%lf %lf %lf ",
		   BottleTime, BottStr1, BottStr2);
	    printf("%u %u %u %lf %lf %lf ",
		   taxonPairDat.numPerTaxa,
		   taxonPairDat.sample[0], taxonPairDat.sample[1],
		   taxonPairDat.tstv[0], taxonPairDat.tstv[1],
		   taxonPairDat.gamma);
	    printf("%u %u %lf %lf %lf ",   
		   taxonPairDat.seqLen, numTauClasses,
		   N1, N2, Nanc);
	    printf("%lf %lf %lf %lf ",
		   taxonPairDat.freqA, taxonPairDat.freqC, taxonPairDat.freqG, 
		   taxonPairDat.freqT);
	    
	    printf("%u\n",  gParam.numTaxaPair);  
	      
	    /* These feed into the system command line (msDQH) within
	            the perl shell msbayes.  Some of these are used directly
		         by msDQH, but some are also passed on to the sumstats
			      programs via the msDQH commabnd line, .... like bp[u],
			      theta, gaussTime, NumPerTax[u], yy, */
	  }
	}

      fprintf(fpTauPsiArray, "%d",numTauClasses);
      if (gParam.numTauClasses > 0) {  /* constrained psi analysis */
for (zzz = 0; zzz < numTauClasses; zzz++)
					fprintf(fpTauPsiArray, "\t%lf",tauArray[zzz]);
					for (zzz = 0; zzz < numTauClasses; zzz++)
					fprintf(fpTauPsiArray, "\t%d",PSIarray[zzz]);
					} 
					fprintf(fpTauPsiArray, "\n");
					}

  fclose(fpTauPsiArray);
  free(tauArray);
  free(PSIarray);
  exit (0);
}
