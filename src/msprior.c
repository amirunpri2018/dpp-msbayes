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
* This program 
*   msprior --reps reps --config file
* reps is number of replications
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
 * Tue Feb 21 2012 Naoki Takebayashi <ntakebayashi@alaska.edu>
 - Revised prior range for Nanc

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

/*
 * Summary of prior distributions, created.
 * [low,high) means uniform distribution with low <= x < high.
 *
 * rec: [0, upperRec * (seqLen-1))
 * mig: [0, upperMig)
 *
 * -- Pop. size (demography) related
 * BottStr1 & 2: [0.01, 1.0) * (N1 or N2)
 *   This assumes the pop size (during the bottleneck) was smaller than current
 * N1 and N2:   (0.01 to 1.99)  (constrained to be N1 + N2 = 2)
 * Nanc  [NancLower/spTheta, gParam.upperAncPopSize * gParam.upperTheta/spTheta)
 *  When there is mut var among loci (THIS IS NOT RIGHT:), 
 * seqLen * [NancLower/locTheta, gParam.upperAncPopSize * gParam.upperTheta/locTheta)
 *     NancLower = max(0.00001 * gParam.lowerTheta, 4*10^(-11))
 *
 * spTheta:  [lowerTheta, upperTheta)       (theta per site)
 *              spTheta <= 0 is excluded even if lowerTheta = 0 or negative.
 * locTheta: spTheta * seqLen * NScaler * mutScaler * mutVar (theta per gene)
 *
 * -- time related
 * tauequalizer = upperTheta / (2 * spTheta * NScaler)
 *
 * # close to 0 means that pop. hasn't started to expand until recently.
 * Bottletime [0.000001, 1.0) * 0.95 * tauequalizer * gaussTime
 *      = scaledGaussTime * [9.5e-7, 0.95) * upperTheta / (2 * spTheta * NScaler)
 *     Some weird lower bound is used.
 *
 * scaledGaussTime: [0.0, upperTau * upperTheta / (2 * spTheta * NScaler))
 *     from taxonTauArray
 *    Then weird lower bound
 */

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>  /* for DBL_EPSILON */

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

int UniqueDouble (double *input, double *output, int inputSize, double smallVal);

int
main (int argc, char *argv[])
{
  double N1, N2, Nanc, NancLower, *uniqTauArray = NULL, *taxonTauArray = NULL,
         *descendant1ThetaArray = NULL, *descendant2ThetaArray = NULL,
         *ancestralThetaArray = NULL, spTheta, thetaMean, descendant1Theta,
         descendant2Theta, tauequalizer, gaussTime = 0.0, mig, rec, BottStr1,
         BottStr2, BottleTime, concentrationParameter;
  double *recTbl;
  int tauClass, *PSIarray = NULL, i, j;
  i_array_2d * divModels;
  i_array * divMod;
  i_array * divIndices;
  int divModelIndex = 0; 
  unsigned int numTauClasses = -1, u, locus, taxonID, zzz;
  unsigned long randSeed;
  unsigned long long rep;
  extern const gsl_rng *gBaseRand;
  int comp_nums (const void *, const void *);

  int b_constrain = 0;
  int *subParamConstrainConfig = NULL;

#ifndef HOMOGENEOUS_MUT
  double *mutScalerTbl;
#endif

  /* set up gParam and gMutParam, as well as gConParam if constrain */
  LoadConfiguration (argc, argv);

  /* set the lower Nanc */
  /* NancLower = 0.00001 * gParam.lowerTheta; */
  /* if (NancLower < 0.00000000004) { /1* 4 * (mu=10^(-11)) * (Ne=1) *1/ */
  /*   NancLower = 0.00000000004; */
  /* } */

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

      for (i = 0; i < strlen (gParam.subParamConstrain); i++)
	{
	  char a = (gParam.subParamConstrain)[i];

	  if (a == '1')
	    subParamConstrainConfig[i] = 1;
	  else if (a == '0')
	    subParamConstrainConfig[i] = 0;
	  else {
	    fprintf(stderr, "ERROR: subParamConstrain string in the config file"
		    "should be either 0 or 1\n");
	    exit (EXIT_FAILURE);
	  }
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

  /* print out all of the parameters */
  if(gParam.printConf) {
    PrintParam(stdout);
    exit (0);
  }

  /* set up arrays */
  /* Sizes are set to the number of taxon pairs (Max number of tau's) */
  if ((b_constrain == 1) && (subParamConstrainConfig[0] == 1)) {
    uniqTauArray = calloc (gParam.numTaxonLocusPairs, sizeof (double));
    PSIarray = calloc (gParam.numTaxonLocusPairs, sizeof (int));
    taxonTauArray = calloc(gParam.numTaxonLocusPairs, sizeof (double));
  } else {
    uniqTauArray = calloc (gParam.numTaxonPairs, sizeof (double));
    PSIarray = calloc (gParam.numTaxonPairs, sizeof (int));
    taxonTauArray = calloc(gParam.numTaxonPairs, sizeof (double));
  }
  descendant1ThetaArray = calloc (gParam.numTaxonPairs, sizeof (double));
  descendant2ThetaArray = calloc (gParam.numTaxonPairs, sizeof (double));
  ancestralThetaArray = calloc (gParam.numTaxonPairs, sizeof (double));
  if ((gParam.concentrationShape > 0) && (gParam.concentrationScale > 0))
  {
      divIndices = init_i_array(gParam.numTaxonPairs);
  }
  else if ((gParam.concentrationShape > -1.0) &&
          (gParam.concentrationScale > -1.0))
  {
      divModels = generate_int_partitions(gParam.numTaxonPairs);
  }
  recTbl = calloc (gParam.numLoci, sizeof (double));

  if (uniqTauArray == NULL || PSIarray == NULL || recTbl == NULL ||
          taxonTauArray == NULL || descendant1ThetaArray == NULL ||
          descendant2ThetaArray == NULL || ancestralThetaArray == NULL)
    {
      fprintf (stderr, "ERROR: Not enough memory for uniqTauArray, PSIarray, or recTbl\n");
      exit (EXIT_FAILURE);
    }

  /* deal with num tau classes */
  if (b_constrain == 0 || subParamConstrainConfig[0] != 1)
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
    }  /* when tau is constrained numTauClasses are set later */

  /* deal with the case when tau is constrained */
  if ((b_constrain == 1) && (subParamConstrainConfig[0] == 1)) {
    int jj, kk;
    double *tempTauArray;
    if ((tempTauArray = calloc(gParam.numTaxonLocusPairs, sizeof(double))) 
	== NULL) {
      fprintf (stderr, "ERROR: Not enough memory for tempTauArray\n");
      exit (EXIT_FAILURE);
    }
    for (jj = 0; jj < gParam.numTaxonLocusPairs; jj++) {
      tempTauArray[jj] = (gConParam.conData[jj]).conTau;
    }
    numTauClasses = UniqueDouble(tempTauArray, uniqTauArray, 
			   gParam.numTaxonLocusPairs, DBL_EPSILON);
    
    if (gParam.numTauClasses != numTauClasses) {
      fprintf (stderr, "WARN: tau's are constrained and found %u different "
	       "classes in the constrain table. But numTauClasses = %u was set."
	       " Using the value found in the constrain table.\n", numTauClasses,
	       gParam.numTauClasses);
      gParam.numTauClasses = numTauClasses;
    } 
    
    /* count tau's to create PSIarray */
    for (jj = 0; jj < gParam.numTaxonLocusPairs; jj++) {
      PSIarray[jj] = 0;
    }
    for (jj = 0; jj < gParam.numTaxonLocusPairs; jj++) {
      for (kk = 0; kk < numTauClasses; kk++) {
	/* there shouldn't be fabs() below */
	if (tempTauArray[jj] - uniqTauArray[kk] < DBL_EPSILON) {
	  PSIarray[kk]++;
	  break;
	}
      }
    }
    free (tempTauArray);
  }

#ifndef HOMOGENEOUS_MUT
  if ((mutScalerTbl = calloc(gParam.numLoci, sizeof(double))) == NULL) {
    fprintf (stderr, "ERROR: Not enough memory for mutScalerTbl\n");
    exit(EXIT_FAILURE);
  }
#endif

  thetaMean = gParam.thetaShape * gParam.thetaScale;

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
            if ((gParam.concentrationShape > 0) &&
                    (gParam.concentrationScale > 0))
            {
                concentrationParameter = gsl_ran_gamma(gBaseRand,
                        gParam.concentrationShape, gParam.concentrationScale);
                numTauClasses = dirichlet_process_draw(gBaseRand,
                        gParam.numTaxonPairs, concentrationParameter, divIndices);
            }
            else if ((gParam.concentrationShape > -1.0) &&
                    (gParam.concentrationScale > -1.0))
            {
                divModelIndex = gsl_rng_uniform_int(gBaseRand,
                        divModels->length);
                divMod = get_i_array_2d(divModels, divModelIndex);
                PSIarray = divMod->a;
                numTauClasses = divMod->length;
            }
            else
            {
            numTauClasses = 1 + gsl_rng_uniform_int (gBaseRand,
                    gParam.numTaxonPairs);
            }
        }
      
      /* create the recombination rate table for each gene */
      rec = 0.0;
      if ((gParam.recombinationShape > 0) && (gParam.recombinationScale > 0))
      {
          rec = gsl_ran_gamma(gBaseRand, gParam.recombinationShape,
                  gParam.recombinationScale);
      }
      /* rec = gsl_ran_flat (gBaseRand, 0.0, gParam.upperRec); */
      for (u=0; u < gParam.numLoci; u++)
	{
	  /* all loci shares same recombination rate */
	  recTbl[u] = rec;
	  /* each locus has different recomb. rate 
	     recTbl[u] = gsl_ran_flat (gBaseRand, 0.0, gParam.upperRec);
	  */
	}
      
#ifndef HOMOGENEOUS_MUT
      /* create regional heterogeneity in the mutation rate */
      if (gParam.numLoci > 1) {
	double shape, scale;
	
	/* arbitrary sample the shape parameter from uniform dist'n */
	shape = gsl_ran_flat(gBaseRand, 1.0, 20);
	/* shape = 1 is exponential with lambda=1, 
	   larger shape -> normal dist'n with smaller var */
	scale = 1/shape; /* E[x] = 1, Var[x] = shape * scale^2 = 1/shape */
	
	/* use gamma */
	for (u=0; u < gParam.numLoci; u++) {
	  mutScalerTbl[u] = gsl_ran_gamma(gBaseRand, shape, scale);
	}
      } else {
	mutScalerTbl[0] = 1.0;
      }
#endif

      // Randomly generate TauArray only when NOT constrain
    if ((b_constrain == 0) || (subParamConstrainConfig[0] != 1))
    {
        int counter;
	    /* sample tau's from uniform prior dist'n */
        for (u = 0; u < numTauClasses; u++)
        {
        // JRO - modified - 11/17/2011
        // uniqTauArray[u] = gsl_ran_flat (gBaseRand, 0.0, gParam.upperTau);
	    /* uniqTauArray[u] = gsl_ran_flat (gBaseRand, gParam.lowerTau, */
	    /*                                 gParam.upperTau); */
            uniqTauArray[u] = gsl_ran_gamma(gBaseRand, gParam.tauShape,
                    gParam.tauScale);
        }

        if ((gParam.concentrationShape > 0) && (gParam.concentrationScale > 0))
        {
            /* memset(PSIarray, 0, sizeof(PSIarray)); */
            for (i = 0; i < gParam.numTaxonPairs; i++)
            {
                PSIarray[i] = 0;
            }
            for (i = 0; i < gParam.numTaxonPairs; i++)
            {
                tauClass = get_i_array(divIndices, i);
                taxonTauArray[i] = uniqTauArray[tauClass];
                PSIarray[tauClass] += 1;
            }
        }
        else if ((gParam.concentrationShape > -1.0) &&
                    (gParam.concentrationScale > -1.0))
        {
            counter = 0;
            for (i = 0; i < divMod->length; i++)
            {
                for (j = 0; j < get_i_array(divMod, i); j++)
                {
                    taxonTauArray[counter] = uniqTauArray[i];
                    counter += 1;
                }
            }
            gsl_ran_shuffle(gBaseRand, taxonTauArray, 
                    gParam.numTaxonPairs, sizeof (double));
        }
        else 
        {
            qsort(uniqTauArray, numTauClasses, sizeof(double),comp_nums);
            for (counter = 0; counter < numTauClasses; counter++) 
            {
                taxonTauArray[counter] = uniqTauArray[counter];
                PSIarray[counter] = 1;
            }

            for (counter = numTauClasses; counter < gParam.numTaxonPairs; counter++)
            {
                tauClass = gsl_rng_uniform_int(gBaseRand, numTauClasses);
                taxonTauArray[counter] = uniqTauArray[tauClass];
                PSIarray[tauClass] = PSIarray[tauClass] + 1;
            }

            /* randomly shuflling the order of taxonTauArray */
            gsl_ran_shuffle(gBaseRand, taxonTauArray, 
                    gParam.numTaxonPairs, sizeof (double));
        }
    }
      
      for (taxonID = 0; taxonID < gParam.numTaxonPairs; taxonID++)
	{
	  //Check upperAncPopSize before doing anything
	  /* ancestral population size prior */
	  /* if (gParam.upperAncPopSize < gParam.lowerTheta) */
	  /*   { */
	  /*     fprintf (stderr, */
		       /* "The upper bound (%lf * %lf) of ancestral pop. size is " */
		       /* "smaller than the lower bound (%lf)\n", */
		       /* gParam.upperAncPopSize, gParam.upperTheta, gParam.lowerTheta); */
	  /*     exit (EXIT_FAILURE); */
	  /*   } */

	  constrainedParameter conTaxonPairDat;

	  /* Population sizes during the bottleneck after the divergence of 2 
	     pops. This is same as the population sizes, immediately after the 
	     divergence/separation of the 2 pops. These are relative sizes. */
	  /* BottStr1 = gsl_ran_flat (gBaseRand, 0.01, 1.0); */
	  /* BottStr2 = gsl_ran_flat (gBaseRand, 0.01, 1.0); */
      BottStr1 = 1.0;
      BottStr2 = 1.0;
      if ((gParam.bottleProportionShapeA > 0) && (gParam.bottleProportionShapeB > 0))
      {
          BottStr1 = gsl_ran_beta(gBaseRand, gParam.bottleProportionShapeA,
                  gParam.bottleProportionShapeB);
          BottStr2 = BottStr1;
          if (gParam.bottleProportionShared == 0) {
              BottStr2 = gsl_ran_beta(gBaseRand, gParam.bottleProportionShapeA,
                      gParam.bottleProportionShapeB);
          }
      }

	  /* After the populations diverge, they experience pop. bottleneck.
	     Then the population size exponentially grows until current size.
	     BottleTime indicate the time when population started to grow.  
	     BottleTime of 1 means, populations start to expand immediately
	     after divergence. Closer to 0 means, populations hasn't started
	     to expand until very recently.  */
	  BottleTime = gsl_ran_flat (gBaseRand, 0.000001, 1.0);

	  /* migration rate prior */
      mig = 0.0;
      if ((gParam.migrationShape > 0) && (gParam.migrationScale > 0))
      {
          mig = gsl_ran_gamma(gBaseRand, gParam.migrationShape,
                  gParam.migrationScale);
      }
	  /* mig = gsl_ran_flat (gBaseRand, 0.0, gParam.upperMig); */

	  /* spTheta prior */

      descendant1Theta = gsl_ran_gamma(gBaseRand, gParam.thetaShape,
              gParam.thetaScale);
      if (gParam.thetaParameters[1] == '1') {
          descendant2Theta = gsl_ran_gamma(gBaseRand, gParam.thetaShape,
                  gParam.thetaScale);
      } else if (gParam.thetaParameters[1] == '0') {
          descendant2Theta = descendant1Theta;
      } else {
          fprintf(stderr, "ERROR: second character of `thetaParameters` "
                "in the config file should be either 0 or 1\n");
	      exit (EXIT_FAILURE);
	  }
      spTheta = (descendant1Theta + descendant2Theta) / 2;
	  /* while ((spTheta = gsl_ran_flat (gBaseRand, gParam.lowerTheta, */
					  /* gParam.upperTheta)) <= 0); */

	  /* The ratio of current population sizes.  The populations
	     exponentially grow to these sizes after bottkleneck is done. */
	  /* both ends excluded for symmetry */
      N1 = descendant1Theta / spTheta;
      N2 = descendant2Theta / spTheta;
	  /* while ((N1 = gsl_ran_flat (gBaseRand, 0.01, 1.99)) == 0.01) */
	  /*   ; */
	  
	  /* N2 = 2.0 - N1; */

	  /* The upper limit of ancestral theta is defined by the product
	     of upper Theta (e.g. 40) and upper AncPopSize (e.g. 0.5) */
	  /* JRO - changing the following hard coded lower limit on ancestral
	     theta to the lower limit specified by user */
	  /* Nanc = gsl_ran_flat (gBaseRand, 0.01,
			       gParam.upperAncPopSize * gParam.upperTheta);*/
	  /* Nanc = gsl_ran_flat (gBaseRand, gParam.lowerTheta, */
			       /* gParam.upperAncPopSize * gParam.upperTheta); */
      if ((gParam.thetaParameters[2] == '2') ||
              ((gParam.thetaParameters[2] == '1') &&
                      (gParam.thetaParameters[1] == '0')))
      {
          if ((gParam.ancestralThetaShape > 0) && (gParam.ancestralThetaScale > 0))
          {
              Nanc = gsl_ran_gamma(gBaseRand, gParam.ancestralThetaShape,
                      gParam.ancestralThetaScale);
          } else
          {
              Nanc = gsl_ran_gamma(gBaseRand, gParam.thetaShape,
                      gParam.thetaScale);
          }
      }
      else if ((gParam.thetaParameters[2] == '1') && 
              (gParam.thetaParameters[1] == '1'))
      {
          Nanc = descendant2Theta;
      }
      else if (gParam.thetaParameters[2] == '0')
      {
          Nanc = descendant1Theta;
      }
      else
      {
          fprintf(stderr, "ERROR: third character of `thetaParameters` "
                "in the config file should be either 0, 1, or 2\n");
	      exit (EXIT_FAILURE);
	  }
      descendant1ThetaArray[taxonID] = descendant1Theta;
      descendant2ThetaArray[taxonID] = descendant2Theta;
      ancestralThetaArray[taxonID] = Nanc;
	  
	  /* pick a tau for every taxon-pair with replacement from the
	     array of X taxon-pairs, where X is a uniform discrete RV
	     from 1 to number of taxon-pairs */
	  if ((b_constrain == 0) || (subParamConstrainConfig[0] != 1))
	    {
	      gaussTime = taxonTauArray[taxonID];
	    }

	  /* use the following if simulating a particular fixed history */
	  /* gaussTime = uniqTauArray[taxonID]; */
	  
	  /* print out the results by going through each locus */
	  for (locus = 0; locus < gParam.numLoci; locus++)
	    {
	      double locTheta, thisNanc, scaledGaussTime, scaledBottleTime;
	      /* check if this locus exist for this taxon pair */
	      /* this table contains 0-offset index for corresponding 
		 taxon:locus mutPara */
	      int mpIndex = gMutParam.locTbl->tbl[taxonID][locus];
	      
	      if(mpIndex<0) { /* this taxon:locus is not in the data */
		continue;
	      }

	      if (b_constrain == 1)
		{  /* If constrained, override with the fixed paras */
		  /* This part is not debugged well 2/14/2008, Naoki */
		  int mpIndex = gMutParam.locTbl->tbl[taxonID][locus];
		  conTaxonPairDat = gConParam.conData[mpIndex];

		  /* tau */
		  /* This allow that tau could differ between loci
		     within a single taxon pair */
		  if (subParamConstrainConfig[0] == 1)
		    gaussTime = conTaxonPairDat.conTau;

		  /** bottleneck priors **/
		  /* severity of bottle neck (how small the pop become) */
		  /* these should be [0,1] */
		  if (subParamConstrainConfig[1] == 1)
		    BottStr1 = conTaxonPairDat.conBottPop1;
		  if (subParamConstrainConfig[2] == 1)
		    BottStr2 = conTaxonPairDat.conBottPop2;
		  
		  /* timing of bottle neck */
		  /* should be [0,1] */
		  if (subParamConstrainConfig[3] == 1)
		    BottleTime = conTaxonPairDat.conBottleTime;
		  
		  /* migration rate prior */
		  if (subParamConstrainConfig[4] == 1)
		    mig = conTaxonPairDat.conMig;
		  
		  /* theta per site */
		  if (subParamConstrainConfig[5] == 1)
		    spTheta = conTaxonPairDat.conTheta;
		  
		  /* population sizes immediately after the separation, and 
		     what it grows to after the bottleneck (today) */
		  /* (0.01, 1.99) */
		  if (subParamConstrainConfig[6] == 1) {
		    N1 = conTaxonPairDat.conN1;
		    N2 = 2.0 - N1;
		  }
		  
		  /* The upper limit of ancestral theta is defined by the 
		     product of upper Theta (e.g. 40) and upper 
		     AncPopSize (e.g. 0.5), then converted to relative size 
		     to spTheta */
		  if (subParamConstrainConfig[7] == 1)
		    /* Nanc = conTaxonPairDat.conNanc * gParam.upperTheta; */
            /* Constrained values of ancestral theta will now have to be
             * specified in absolute terms */
		    Nanc = conTaxonPairDat.conNanc;
		  
		  /* recombination rate per neighboring site */
		  if (subParamConstrainConfig[8] == 1)
		    recTbl[locus] = conTaxonPairDat.conRec;
		}  /* end of constrai */

	      /* access sample sizes, mutational model for this taxon:locus */
	      mutParameter taxonPairDat;
	      taxonPairDat = gMutParam.data[mpIndex];
	      
	      /* scale the theta for each locus */
	      /* Note that species wide theta (represents pop size) is 
	         4 Ne mu with mu per site, not per gene.
		 Assumes mu is constant.  This may be a problem with
	         mitochondoria */
	      locTheta = spTheta * taxonPairDat.seqLen * 
		taxonPairDat.NScaler * taxonPairDat.mutScaler;
#ifndef HOMOGENEOUS_MUT
	      locTheta *=  mutScalerTbl[locus];
#endif

	      /* thisNanc is basically a random deviate from a uniform dist'n:
		 [gParam.lowerTheta / spTheta, 
		   gParam.upperAncPopSize * gParam.upperTheta/spTheta) 
		 For example, if upperTheta = 10 & upperAncPopSize = 0.5,
		 upperAncTheta become 10 * 0.5 = 5.
		 msDQH specify the past population sizes in terms of the 
		 ratio of N_anc / N_theta, so the following division
		 by locTheta is required.
	      */
	      /* thisNanc = Nanc * taxonPairDat.seqLen / locTheta; */
	      thisNanc = Nanc / spTheta; /* this can be done outside of locus loop */

	      /* this scaling is done inside of locus loop to accomodate 
		 the gamma dist'n of mut rate for each locus */

	      /* tauequalizer = gParam.upperTheta / */ 
		/* 2 / (spTheta * taxonPairDat.NScaler); */
          tauequalizer = thetaMean / (spTheta * taxonPairDat.NScaler);
	      /* WORK, CONFIRM THIS. Naoki Nov 2, 2009.  IT USED TO BE
		 tauequalizer = gParam.upperTheta * taxonPairDat.seqLen / 
		 2 / locTheta;

	      */

	      /* Division by 2 is coming from N1 + N2 = 2.
		 We are considering that N_0 in theta_0 (=4 N_0 mu) specified 
		 for -t option (we use -t locTheta) of msDQH is equal to 
		 (N1+N2)/2 */

	      scaledGaussTime = gaussTime * tauequalizer;
	      /* 1 unit of tau (gaussTime) = 2 N_max (N_max is the 
		 N assumed in upperTheta) */
	      /* I think we should get rid of /2 from tauequalizer */

	      /* The following if is a little weird */
	      if (scaledGaussTime < 0.0001) {
		scaledGaussTime  = 0.0001;
		scaledBottleTime = 0.00005;
	      } else {
		scaledBottleTime = BottleTime * 0.95 * scaledGaussTime;
	      }
	      
	      if (debug_level)
		fprintf (stderr, 
			 "DEBUG: scaled BottleTime:%lf\tgaussTime:%lf\n",
			 scaledBottleTime, scaledGaussTime);

	      /* We can send some extra info to msbayes.pl here */
	      printf ("%u %u %u ", lociTaxonPairIDcntr, taxonID+1, locus+1);
	      lociTaxonPairIDcntr ++; /* seriral id: 1 to # taxon:locus pairs */
	      printf ("%.11lf %.11lf %.11lf %.11lf ",
		      locTheta, scaledGaussTime, mig, 
		      recTbl[locus] * (taxonPairDat.seqLen - 1));
	      printf ("%.11lf %.11lf %.11lf ", scaledBottleTime, 
		      BottStr1 * N1, BottStr2 * N2);
	      printf ("%u %u %u %lf %lf %lf ",
		      taxonPairDat.numPerTaxa,
		      taxonPairDat.sample[0], taxonPairDat.sample[1],
		      taxonPairDat.tstv[0], taxonPairDat.tstv[1],
		      taxonPairDat.gamma);
	      printf ("%u %.11lf %.11lf %.17lf ",
		      taxonPairDat.seqLen, N1, N2, thisNanc);
	      printf ("%lf %lf %lf %lf\n",
		      taxonPairDat.freqA, taxonPairDat.freqC,
		      taxonPairDat.freqG, taxonPairDat.freqT);

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
	printf (",%lf", uniqTauArray[zzz]);
      printf(" psiTbl:");
      for (zzz = 0; zzz < numTauClasses; zzz++)
	printf (",%d", PSIarray[zzz]);
      printf(" d1ThetaTbl:");
      for (zzz = 0; zzz < gParam.numTaxonPairs; zzz++)
          printf (",%lf", descendant1ThetaArray[zzz]);
      printf(" d2ThetaTbl:");
      for (zzz = 0; zzz < gParam.numTaxonPairs; zzz++)
          printf (",%lf", descendant2ThetaArray[zzz]);
      printf(" aThetaTbl:");
      for (zzz = 0; zzz < gParam.numTaxonPairs; zzz++)
          printf (",%lf", ancestralThetaArray[zzz]);
      printf("\n");

    }

  free (uniqTauArray);
  free (taxonTauArray);
  free (PSIarray);
  free (descendant1ThetaArray);
  free (descendant2ThetaArray);
  free (ancestralThetaArray);
  free (recTbl);
  free (subParamConstrainConfig);
  if ((gParam.concentrationShape > 0) && (gParam.concentrationScale > 0))
  {
      free_i_array(divIndices);
  }
  else if ((gParam.concentrationShape > -1.0) &&
          (gParam.concentrationScale > -1.0))
  {
      free_i_array_2d(divModels);
  }
  exit (0);
}


/* 
 * Arguments:
 *   input: array with inputSize elements (will not be modified)
 *   output: result will be returned in this array.
 *           the memory for inputSize elements should be allocated
 *   inputSize: number of elements in input array
 *   smallVal: used to decide whether two double number differs or not.
 *      DBL_EPSILON defined in float.h may be appropriate.
 * 
 * Extract the unique elements and the unique array will be returned.
 * The result array is sorted from the smallest to largest.
 *
 * If the difference between two numbers is < smallVal, the two numbers 
 * are considered as identical.
 *
 * Example: input = (0.1, 0.4, 0.2, 0.3, 0.51) and smallVal = 0.1.
 * It will return (0.1, 0.3, 0.51).
 *
 * Explanation:
 * It sorts the array at first, and then works from the smallest number.
 * 0.2 is removed (=0.1).  Since 0.2 is removed, 0.3 is
 * not equal to 0.1 (so it will be retained).
 *
 * The number of elements in the unique array is returned.  Elements after
 * the unique values aren't initialized (contains some numbers from input).
 */

int
UniqueDouble (double *input, double *output, int inputSize, double smallVal) {
  int checkI;
  int lastIndex = 0;

  if (input == NULL) {
    fprintf (stderr, "ERROR, input array is empty in UniqueDouble\n");
    return -1;
  }

  if (inputSize == 0) {
    return 0;
  }

  if (input != output) {
    memcpy (output, input, inputSize * sizeof (double));
  }
  qsort(output, inputSize, sizeof(double), comp_nums);

  for (checkI = 1; checkI < inputSize; checkI++) {
    if (output[checkI] - output[lastIndex] >= smallVal) {
      /* sufficiently different */
      output[lastIndex+1] = output[checkI];
      lastIndex++;
    }
  }

  return (lastIndex + 1);
}
