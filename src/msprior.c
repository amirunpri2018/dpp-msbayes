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
 * Fri July 26 2006 Mike Hickerson
 - Implement fixed number of tau classes.

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
#include <string.h>

#include <gsl/gsl_rng.h>	/* for base rand num gen's */
#include <gsl/gsl_randist.h>	/* for gsl_ran_gamma */
#include "msprior.h"
#include "setup.h"

/* This has to be Global */
const gsl_rng *gBaseRand;	/* global rand number generator */

/* 
* global variables which stores the parameters (settings for the upper
* limist of prior dist'n etc
*/
runParameters gParam;		/* stores upper limits of prior dist'n etc */
mutParameterArray gMutParam;	/* stores mut model & # samples for 
				   each taxon pair */
constrainedParameterArray gConParam;	/* store constrained sub-parameters for
					   each taxon pair */


int
comp_nums (const void *doubleNum1, const void *doubleNum2)
{
  const double *num1 = doubleNum1;
  const double *num2 = doubleNum2;
  if (*num1 < *num2)
    return -1;
  else if (*num1 > *num2)
    return 1;
  else
    return 0;
}

int
main (int argc, char *argv[])
{
  double N1, N2, Nanc, *tauArray = NULL, spTheta, tauequalizer, gaussTime = 0.0,
    mig, rec, BottStr1, BottStr2, BottleTime;
  double *recTbl;
  int tauClass, *PSIarray = NULL;
  unsigned int numTauClasses = -1, u, locus, taxonID, zzz, c;
  unsigned long randSeed;
  unsigned long long rep;
  extern const gsl_rng *gBaseRand;
  int comp_nums (const void *, const void *);
  FILE *fpTauPsiArray;

  int b_constrain = 0;
  int *subParamConstrainConfig = NULL;

  /* set up gParam and gMutParam, as well as gConParam if constrain */
  LoadConfiguration (argc, argv);

  /* set b_constrain to 1 if constrain */
  if (gParam.constrain > 0)
    {
      //initialize constrain indicator
      b_constrain = 1;


      //initialize subParamConstrainConfig array
      subParamConstrainConfig = calloc (NUMBER_OF_CONPARAM, sizeof (int));
      if (subParamConstrainConfig == NULL)
	{
	  fprintf (stderr,
		   "ERROR: Not enough memory for subParamConstrainConfig\n");
	  exit (EXIT_FAILURE);
	}

      int i = 0;

      for (i = 0; i < strlen (gParam.subParamConstrain); i++)
	{
	  char a = (gParam.subParamConstrain)[i];

	  if (a == '1')
	    subParamConstrainConfig[i] = 1;
	  else if (a == '0')
	    subParamConstrainConfig[i] = 0;
	}

    }

  /* for initiating the gsl random number generator */
  /* initialize PRNG */
  srand (gParam.prngSeed);	/* Better way of seeding here ? */
  randSeed = rand ();
  if (debug_level > 0)
    randSeed = 1;

  gBaseRand = gsl_rng_alloc (gsl_rng_mt19937);	/* set the base PRNG to
						   Mersenne Twister */
  gsl_rng_set (gBaseRand, randSeed);	/* seed the PRNG */

  if (b_constrain == 0)
    {
      /* fixed numTauClasses configuration */
      if (gParam.numTauClasses != 0)
	{
	  if (gParam.numTauClasses > gParam.numTaxonPairs)
	    {
	      fprintf (stderr, "WARN: numTauClasses (%u) is larger than "
		       "numTaxonPairs (%u). Setting numTauClasses to %u",
		       gParam.numTauClasses, gParam.numTaxonPairs,
		       gParam.numTaxonPairs);
	      gParam.numTauClasses = gParam.numTaxonPairs;
	    }
	  numTauClasses = gParam.numTauClasses;
	}
    }
  else
    {
      numTauClasses = gParam.numTauClasses;
    }

  /* print out all of the parameters */
  if(gParam.printConf) {
    PrintParam(stdout);
    exit (0);
  }

  /* Sizes are set to the number of taxon pairs (Max number of tau) */
  tauArray = calloc (gParam.numTaxonPairs, sizeof (double));
  PSIarray = calloc (gParam.numTaxonPairs, sizeof (int));

  recTbl = calloc (gParam.numLoci, sizeof (double));
  if (tauArray == NULL || PSIarray == NULL || recTbl == NULL)
    {
      fprintf (stderr, "ERROR: Not enough memory for tauArray, PSIarray, or recTbl\n");
      exit (EXIT_FAILURE);
    }
  
  /* open the stream to output tauArray and PSIarray */
  if ((fpTauPsiArray = fopen (gParam.priorOutFile, "w")) == NULL)
    {
      fprintf (stderr, "Cannot open the file.\n");
      exit (1);
    }

  /* print out the column headers to the tauArray and PSIarray  */
  /* PRI.numTauClass PRI.Tau.1  PRI.Tau.2 PRI.Tau.3 ... PRI.Psi.1 PRI.Psi.2 PRI.Psi.3 ... */
  fprintf (fpTauPsiArray, "PRI.numTauClass");
  if (gParam.numTauClasses > 0)
    {				/* constrained psi analysis */
      for (zzz = 0; zzz < numTauClasses; zzz++)
	fprintf (fpTauPsiArray, "\tPRI.Tau.%d", zzz + 1);
      for (zzz = 0; zzz < numTauClasses; zzz++)
	fprintf (fpTauPsiArray, "\tPRI.Psi.%d", zzz + 1);
    }
  fprintf (fpTauPsiArray, "\n");


  /* Beginning of the main loop */
  for (rep = 0; rep < gParam.reps; rep++)
    {
      int lociTaxonPairIDcntr = 1;
      /*
       * Each taxon pair was separated at a time tau in the past.  Of
       * all pairs, some of them may have been separated at the same
       * time.  numTauClasses is the number of classes with different
       * divergence time.
       *
       * If gParam.numTauClasses is not set, we are sampling
       * numTauClasses from a uniform prior dist'n.
       */
      if (gParam.numTauClasses == 0)
	{			/* numTauClasses is NOT fixed */
	  numTauClasses =
	    1 + gsl_rng_uniform_int (gBaseRand, gParam.numTaxonPairs);
	}
      
      for (c = 0; c < numTauClasses; c++)
	PSIarray[c] = 0;	/* Reset the PSIarray counters */

      int psiIndex = 0;
      // Randomly generate TauArray only when NOT constrain
      if ((b_constrain == 0) || (subParamConstrainConfig[0] != 1))
	{
	  /* sample tau's from uniform prior dist'n */
	  for (u = 0; u < numTauClasses; u++)
	    {
	      tauArray[u] = gsl_ran_flat (gBaseRand, 0.0, gParam.upperTau);

	      if (debug_level)
		{
		  fprintf (stderr, "DEBUG:%u of %u categories:\t%lf\n",
			   u, numTauClasses, tauArray[u]);
		}
	    }

	  /* create the recombination rate table for each gene */
	  rec = gsl_ran_flat (gBaseRand, 0.0, gParam.upperRec);
	  for (u=0; u <gParam.numLoci; u++)
	    {
	      /* all loci shares same recombination rate */
	      recTbl[u] = rec;
	      /* each locus has different recomb. rate 
	      recTbl[u] = gsl_ran_flat (gBaseRand, 0.0, gParam.upperRec);
	      */
	    }
	}
      else if ((b_constrain == 1) && (subParamConstrainConfig[0] == 1))
	{

	  for (u = 0; u < numTauClasses; u++)
	    tauArray[u] = (gConParam.conData[u]).conTau;

	  double *conTauArray = NULL;

	  conTauArray = calloc (gParam.numTaxonPairs, sizeof (double));
	  if (conTauArray == NULL)
	    {
	      fprintf (stderr, "ERROR: Not enough memory for conTauArray\n");
	      exit (EXIT_FAILURE);
	    }

	  for (u = 0; u < gParam.numTaxonPairs; u++)
	    {
	      conTauArray[u] = gConParam.conData[u].conTau;
	    }

	  qsort (conTauArray, (gParam.numTaxonPairs), sizeof (double),
		 comp_nums);


	  //initialize PSIarray
	  double tempTau = conTauArray[0];

	  PSIarray[psiIndex] = 1;
	  for (u = 1; u < gParam.numTaxonPairs; u++)
	    {

	      if (conTauArray[u] == tempTau)
		PSIarray[psiIndex]++;
	      else
		{
		  tempTau = conTauArray[u];
		  PSIarray[++psiIndex]++;
		}
	    }

	  free (conTauArray);
	}

      //if(psiIndex != (numTauClasses-1) )
      // fprintf(stderr,"numberPsiArray:%d and numTauClasses:%d have problems\n", psiIndex, numTauClasses );

      int tauPsiIndex = 0;
      int tauCounter = 1;
      for (taxonID = 0; taxonID < gParam.numTaxonPairs; taxonID++)
	{
	  //Check upperAncPopSize before doing anything
	  /* ancestral population size prior */
	  if (gParam.upperAncPopSize < 0.01)
	    {
	      fprintf (stderr,
		       "The upper bound (%lf) of ancestral pop. size is "
		       "smaller than the lower bound (0.01)\n",
		       gParam.upperAncPopSize);
	      exit (EXIT_FAILURE);
	    }


	  constrainedParameter conTaxonPairDat;

	  /* Population sizes during the bottleneck after the divergence of 2 pops.
	     This is same as the population sizes, immediately after the 
	     divergence/separation of the 2 pops. These are relative sizes. */
	  BottStr1 = gsl_ran_flat (gBaseRand, 0.01, 1.0);
	  BottStr2 = gsl_ran_flat (gBaseRand, 0.01, 1.0);

	  /* timing of bottle neck. This is the time when
	     After the populations diverge, they experience pop. bottleneck.  Then
	     the population size exponentially grow until current.
	     BottleTime indicate the time when population start to grow.  
	     BottleTime of 1 means, populations start to expand immediately after
	     divergence. Closer to 0 means, populations hasn't started to expand
	     until very recently.  */
	  BottleTime = gsl_ran_flat (gBaseRand, 0.000001, 1.0);

	  /* migration rate prior */
	  mig = gsl_ran_flat (gBaseRand, 0.0, gParam.upperMig);
	  /* spTheta prior */
	  spTheta = gsl_ran_flat (gBaseRand, gParam.lowerTheta,
				gParam.upperTheta);

	  /* The ratio of current population sizes.  The populations exponentially
	     grow to these sizes after bottkleneck is done. */
	  N1 = gsl_ran_flat (gBaseRand, 0.01, 1.99);
	  /* WORK: Mike is the upper limit of 1.99 ok, even for nuclear gene? NT Aug 26, 2008*/
	  /*
	     Nmax=((gParam.upperAncPopSize*gParam.upperTheta)*gParam.upperTheta)
	     /spTheta; 
	     Nanc = gsl_ran_flat (gBaseRand, 0.01, Nmax);
	   */

	  /* The upper limit of ancestral theta is defined by the product
	     of upper Theta (e.g. 40) and upper AncPopSize (e.g. 0.5) */
	  Nanc = gsl_ran_flat (gBaseRand, 0.01,
			       gParam.upperAncPopSize * gParam.upperTheta);
	  
	  /* sample sizes, constrain model for taxonID-th taxon-pair */
	  //conTaxonPairDat = gConParam.conData[taxonID];

	  if (b_constrain == 1)
	    {
	      int gInd;
	      conTaxonPairDat = gConParam.conData[taxonID];

	      /** bottleneck priors **/
	      /* severity of bottle neck (how small the population become) */
	      if (subParamConstrainConfig[1] == 1)
		BottStr1 = conTaxonPairDat.conBottPop1;
	      if (subParamConstrainConfig[2] == 1)
		BottStr2 = conTaxonPairDat.conBottPop2;
	      /* timing of bottle neck */
	      if (subParamConstrainConfig[3] == 1)
		BottleTime = conTaxonPairDat.conBottleTime;

	      /* migration rate prior */
	      if (subParamConstrainConfig[4] == 1)
		mig = conTaxonPairDat.conMig;

	      /* theta prior */
	      if (subParamConstrainConfig[5] == 1)
		spTheta = conTaxonPairDat.conTheta;

	      /* population sizes immediately after the separation, and 
	         what it grows to after the bottleneck (today) */
	      if (subParamConstrainConfig[6] == 1)
		N1 = conTaxonPairDat.conN1;

	      /* The upper limit of ancestral theta is defined by the product
	         of upper Theta (e.g. 40) and upper AncPopSize (e.g. 0.5) */
	      if (subParamConstrainConfig[7] == 1)
		Nanc = conTaxonPairDat.conNanc;

	      /* recombination rate */
	      if (subParamConstrainConfig[8] == 1) {
		rec = conTaxonPairDat.conRec;
		/* all loci have the same recombination rate */
		/* check this with Wen and Mike */
		for (gInd = 0; gInd < gParam.numLoci; gInd++) {
		  recTbl[gInd] = rec;
		}
	      }
	    }

	  N2 = 2.0 - N1;
	  /* WORK: Mike, is 2.0 ok, even for nuclear gene? NT Aug 26, 2008*/

	  Nanc = Nanc / spTheta; /* get the ratio of theta_anc / spTheta_cur 
				    This ratio is required for msDQH */

	  tauequalizer = gParam.upperTheta / 2 / spTheta;
	  /* WORK: Mike, is "2" ok for nuclear gene? NT Aug 26, 2008 */

	  /* pick a tau for every taxon-pair with replacement from the
	     array of X taxon-pairs, where X is a uniform discrete RV
	     from 1 to number of taxon-pairs */

	  if ((b_constrain == 0) || (subParamConstrainConfig[0] != 1))
	    {
	      // 1-31-2007 it will use new way of picking tau
	      /* WORK: this method has a bias in multi locus case NT*/
	      if (taxonID < numTauClasses)
		tauClass = taxonID;
	      else
		tauClass = gsl_rng_uniform_int (gBaseRand, numTauClasses);

	      gaussTime = tauArray[tauClass];
	      PSIarray[tauClass] = PSIarray[tauClass] + 1;

	      //printf("picking index: %d, gaussTime: %lf  ", tauClass, tauArray[tauClass]);
	    }
	  else if ((b_constrain == 1) && (subParamConstrainConfig[0] == 1))
	    {
	      if (taxonID < numTauClasses)
		{
		  tauClass = taxonID;
		}
	      else
		{
		  if (tauCounter <= ((PSIarray[tauPsiIndex]) - 1))
		    {
		      tauClass = tauPsiIndex;
		      tauCounter++;
		    }
		  else
		    {
		      tauCounter = 1;
		      int x = 1;
		      for (x = 1; x < numTauClasses; x++)
			if (gConParam.conData[taxonID].conTau == tauArray[x])
			  tauPsiIndex = x;
		      //tauPsiIndex ++;
		      tauClass = tauPsiIndex;
		      tauCounter++;
		    }
		}

	      if (tauClass < numTauClasses)
		gaussTime = tauArray[tauClass];
	      else
		fprintf (stderr, "what happened to tauClass?\n");
	    }
	  //gaussTime = tauArray[taxonID];
	  //fprintf(stderr, "gaussTime")

	  /* use the following if simulating a particular fixed history */
	  /* gaussTime = tauArray[taxonID]; */

	  gaussTime = gaussTime * tauequalizer;

	  /* The following 2 if's are weird */
	  if (gaussTime < 0.0001)
	    gaussTime = 0.0001;

	  BottleTime = BottleTime * 0.95 * gaussTime;

	  if (gaussTime < 0.0001)
	    BottleTime = 0.00005;

	  if (debug_level)
	    fprintf (stderr, "DEBUG: BottleTime:%lf\tgaussTime:%lf\n",
		     BottleTime, gaussTime);

	  /* print out the results */
	  for (locus = 0; locus < gParam.numLoci; locus++)
	    {
	      double locTheta;
	      /* check if this locus exist for this taxon pair */
	      /* this table contains 0-offset index for corresponding 
		 taxon:locus mutPara */
	      int mpIndex = gMutParam.locTbl->tbl[taxonID][locus];
	      
	      if(mpIndex<0) { /* this taxon:locus is not in the data */
		continue;
	      }

	      /* access sample sizes, mutational model for this taxon:locus */
	      mutParameter taxonPairDat;
	      taxonPairDat = gMutParam.data[mpIndex];
	      
	      /* scale the theta for each locus */
	      /* Note that species wide theta (represents pop size) is 
	         4 Ne mu with mu per site, not per gene.
		 Assumes mu is constant.  This may be a problem with
	         mitochondoria */
	      locTheta = spTheta * taxonPairDat.seqLen * taxonPairDat.ploidy/2;
	      
	      /* We can send some extra info to msbayes.pl here */
	      printf ("%u %u %u ", lociTaxonPairIDcntr, taxonID+1, locus+1);
	      lociTaxonPairIDcntr ++; /* seriral id: 1 to # taxon:locus pairs */
	      printf ("%lf %lf %lf %lf ",
		      locTheta, gaussTime, mig, recTbl[locus]);
	      printf ("%lf %lf %lf ", BottleTime, BottStr1, BottStr2);
	      printf ("%u %u %u %lf %lf %lf ",
		      taxonPairDat.numPerTaxa,
		      taxonPairDat.sample[0], taxonPairDat.sample[1],
		      taxonPairDat.tstv[0], taxonPairDat.tstv[1],
		      taxonPairDat.gamma);
	      printf ("%u %lf %lf %lf ",
		      taxonPairDat.seqLen, N1, N2, Nanc);
	      printf ("%lf %lf %lf %lf ",
		      taxonPairDat.freqA, taxonPairDat.freqC,
		      taxonPairDat.freqG, taxonPairDat.freqT);
	      printf ("%u\n", numTauClasses);
	      /* These feed into the system command line (msDQH) within
	         the perl shell msbayes.  Some of these are used directly
	         by msDQH, but some are also passed on to the sumstats
	         programs via the msDQH commabnd line, .... like bp[taxonID],
	         theta, gaussTime, NumPerTax[taxonID], yy, */
	    }
	}

      /* The followings are used to calculate prior, processed in msbayes.pl */
      printf ("# TAU_PSI_TBL setting: %d realizedNumTauClasses: %u tauTbl:", 
	      gParam.numTauClasses, numTauClasses);
      for (zzz = 0; zzz < numTauClasses; zzz++)
	printf (",%lf", tauArray[zzz]);
      printf(" psiTbl:");
      for (zzz = 0; zzz < numTauClasses; zzz++)
	printf (",%d", PSIarray[zzz]);
      printf("\n");

#if 0
      fprintf (fpTauPsiArray, "%d", gParam.numTauClasses);
      if (gParam.numTauClasses > 0)
	{			/* constrained psi analysis */
	  for (zzz = 0; zzz < numTauClasses; zzz++)
	    fprintf (fpTauPsiArray, "\t%lf", tauArray[zzz]);
	  for (zzz = 0; zzz < numTauClasses; zzz++)
	    fprintf (fpTauPsiArray, "\t%d", PSIarray[zzz]);
	}
      fprintf (fpTauPsiArray, "\n");
#endif
    }

  fclose (fpTauPsiArray);
  free (tauArray);
  free (PSIarray);
  free (recTbl);
  free (subParamConstrainConfig);
  exit (0);
}
