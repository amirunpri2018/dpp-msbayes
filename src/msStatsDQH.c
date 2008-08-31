/*
 * msStatsDQH.c
 *
 * Copyright (C) 2006  Richard Hudson, Eli Stahl,
 *                     Michael Hickerson, Naoki Takebayashi
 *
 * This file is a part of sumstatsvector, distributed with msBayes.
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

#include "hashtab.h"
#include "msprior.h"
#include "sumStatsVector.h"
#include <string.h>

/* #include <float.h> */
/* #if defined (__SVR4) && defined (__sun) */
/*   int isinf(double x) { return !finite(x) && x==x; } */
/* #endif */

#define MAX_LEN_COLUMN_NAME 128	/* Used for header. This is the maximum char length
				   of names for each column */

double nucdiv (int, int, char **);
double nucdiv_w (int, int, char **, int, int *, int);
double nucdiv_bw (int, int, char **, int, int *);
double tajddenominator (int, int, double);
double tajddenominator2 (int, int, double);
double thetaW (int, int);
double thetah (int, int, char **);
void FuLi (double *D, double *F, int, int, char **, double pi);
void FrequencyDistrInfSites (int *freqdist, int nsam, int segsites,
			     char **list);
void FrequencyDistrQ (int *freqdist, int nsam, int segsites, char **list);
int segsub (int nsub, int segsites, char **list);
void segsubs (int *segwithin, int segsites, char **list, int npops, int *n);
int multiplepopssampledfrom (int nsam, int npops, int *n);	/*zzz n[] is called config[] in samples.c zzz */
static int SS_comp (const void *, const void *);
#if 0
static int compare_doubles (const void *a, const void *b);
#endif
int frequency (char, int, int, char **);


void shannonIndex (char **list, int *config, double **shannonIndexArray);
int charCount (char *arr);
extern int gPrintHeader;	/* boolean 1 means print header (column names), 0 = no header
				   -H option invoke printing of the header */
static void PrintHeader (char priorNames[][MAX_LEN_COLUMN_NAME],
			 int numPriors,
			 char sumStatNames[][MAX_LEN_COLUMN_NAME],
			 int numSumStats, int numTaxonPairs);

/***** MAKE SURE the following two lines are up-to-date *****/
int numSumStats = 17;
char ssNameVect[][MAX_LEN_COLUMN_NAME] =
  { "pi.b", "pi.w", "pi", "wattTheta", "pi.net", "tajD", "tajD.denom",
  "pi.wPop2", "pi.wPop1", "wattTheta.Pop2", "wattTheta.Pop1",
  "tajD.denomPop2", "tajD.denomPop1", "ShannonsIndex.Between",
  "ShannonsIndex.Net", "ShannonsIndex.Pop1", "ShannonsIndex.Pop2"
};


#if 0				/* commented out since it is not used */
/*
 * used for qsort to compare two doubles.
 * Takes two pointers to double as the arguments.
 *
 * Returns: 1  if a > b
 *          0  if a == b
 *         -1  if a < b
 */

static int
compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}
#endif


static int
SS_comp (const void *p1, const void *p2)
{
  const struct SumStat *sp1 = (struct SumStat *) p1;
  const struct SumStat *sp2 = (struct SumStat *) p2;

  return ((sp1->PI_b) > (sp2->PI_b)) - ((sp1->PI_b) < (sp2->PI_b));
  /*return  ( sp1->PI_b) - ( sp2->PI_b); */
}


/* Print out the available summary stats and Exit */
void
PrintSumStatNames (void)
{
  int i;
  for (i = 0; i < numSumStats; i++)
    {
      printf ("%s\n", ssNameVect[i]);
    }
  exit (0);
}

/*
  Takes a pointer to structure which contains the results of SINGLE run of
  msDQH and calculates the summary statistics.

  Side effect: memory is allocated to store all SumStats,
  Returns: a pointer to the SumStats which contains the results

  nsam: number of samples
  segsites: number of segregating sites
  list: character matrix (A,T,G,orC) containing nsam rows and segsites columns
  nsub: gNadv, default 0, but can be specified by --nadv option
  npops: number of sub-populations
  n: a vector of sub-population sizes, there are npops elements
  theta: 4 Ne mu used for the simulations, it comes from the command line 
         option (-t) to msDQH 
*/
struct SumStat *
CalcSumStats (msOutput *msOut)
{
  int i, nsam, npops, *n, BasePairs, segsites;
  int *freqdist, *segwithin, tW_w_npops;
  double  FuLiD, FuLiF, Fst, Nm, tW_w;
  /* double h, th, ObsvVARD, ObsvVARPi_Net, ObsvEPi_Net, ObsvCV_pi_net_tW; */
  char **list;

  /* struct SumStat SumStat_list[NumTaxa]; */
  struct SumStat * resultSS;

  /* copying frequently used valuse for easy access */
  nsam = msOut->nsam;
  npops = msOut->npops;
  n = msOut->n;
  BasePairs = msOut->BasePairs;
  segsites = msOut->segsites;
  list = msOut->seqDat;

  /* SIDE EFFECT: allocating memory to return results */
  if (!(resultSS = malloc (sizeof(struct SumStat)))) {
    perror("ERROR: No Mem in CalcSumStats\n");
    exit(EXIT_FAILURE);
  }
  
  /* initalize with -1 */
  resultSS->PI_b = resultSS->PI_Net = resultSS->PI_w2 = resultSS->PI_w1 =
    resultSS->PI_w = -1;

  if (! (freqdist = (int *) malloc (nsam * sizeof (int)))) {
	perror("ERROR: No Mem 2 in CalcSumStats\n");
	exit(EXIT_FAILURE);    
  }
  if (msOut->Qbool)
    FrequencyDistrQ (freqdist, nsam, segsites, list);
  else
    FrequencyDistrInfSites (freqdist, nsam, segsites, list);

  /* this is actually already taken care of, so this is redundant */
  if (msOut->isNumSegSitesFixed)
    resultSS->TW = msOut->theta;
  else
    resultSS->TW = thetaW (nsam, segsites);

  resultSS->PI = nucdiv (nsam, segsites, list);

  if (msOut->Fst_bool)
    {
      if (!(segwithin = (int *) malloc (npops * sizeof (int)))) {
	perror("ERROR: No Mem 3 in CalcSumStats\n");
	exit(EXIT_FAILURE);
      }
      
      segsubs (segwithin, segsites, list, npops, n);

      tW_w = 0.;
      tW_w_npops = 0;
      for (i = 0; i < npops; i++)
	if (n[i] > 1)
	  {
	    tW_w += thetaW (n[i], segwithin[i]);
	    tW_w_npops++;
	  }
      tW_w /= tW_w_npops;

      /*yyy ABOVE  Eli 05/15/06 yyy */
      resultSS->TW1 = thetaW (n[0], segwithin[0]);
      resultSS->TW2 = thetaW (n[1], segwithin[1]);
      /* -1 signals Average of pi's within subpop */

      resultSS->PI_w =
	nucdiv_w (nsam, segsites, list, npops, n, -1) / BasePairs;
      resultSS->PI_w1 = nucdiv_w (nsam, segsites, list, npops, n, 0);
      resultSS->PI_w2 = nucdiv_w (nsam, segsites, list, npops, n, 1);

      resultSS->PI_b = nucdiv_bw (nsam, segsites, list, npops, n) / BasePairs;
      Fst = 1. - resultSS->PI_w / resultSS->PI_b;
      resultSS->PI_Net = resultSS->PI_b - resultSS->PI_w;
      if (Fst < 0)
	{
	  Fst = 0.;
	  Nm = -1.;
	}
      else
	{
	  Nm = (1. / Fst - 1.) / 4.;
	}
    }
  /*   th = thetah(nsam, segsites, list) ; */

  /* Tajima's D denominator */
  resultSS->TDD = tajddenominator (nsam, segsites, resultSS->PI);
  resultSS->TDD1 = tajddenominator2 (n[0], segwithin[0], resultSS->PI_w1);
  resultSS->TDD2 = tajddenominator2 (n[1], segwithin[1], resultSS->PI_w2);

  /* Tajima's D */
  resultSS->TD = (resultSS->PI - resultSS->TW) / resultSS->TDD;
  resultSS->TD1 = (resultSS->PI_w1 - resultSS->TW1) / resultSS->TDD1;
  resultSS->TD2 = (resultSS->PI_w2 - resultSS->TW2) / resultSS->TDD2;
  
/*  FuLi(&FuLiD,&FuLiF,nsam,segsites,list,resultSS->PI);*/
/*   h = resultSS->PI-th ; */
#if 0
  /* not used NT */
  if (msOut->nsub > 0)
    nsegsub = segsub (msOut->nsub, segsites, list);
#endif
  resultSS->TW = resultSS->TW / BasePairs;
  resultSS->TW1 = resultSS->TW1 / BasePairs;
  resultSS->TW2 = resultSS->TW2 / BasePairs;

  resultSS->PI = resultSS->PI / BasePairs;
  resultSS->PI_w1 = resultSS->PI_w1 / BasePairs;
  resultSS->PI_w2 = resultSS->PI_w2 / BasePairs;

  if (segsites < 1)
    resultSS->PI_b = resultSS->PI_Net = 
      resultSS->TD = resultSS->TD1 = resultSS->TD2 = 
      resultSS->PI_w = resultSS->PI_w1 = resultSS->PI_w2 = resultSS->PI = 
      resultSS->TW = resultSS->TW1 = resultSS->TW2 =
      resultSS->TDD = resultSS->TDD1 = resultSS->TDD2 = 
      FuLiD = FuLiF = Fst = 0;
  
  if (segwithin[1] < 1)
    resultSS->TD2 = 0;
  if (segwithin[0] < 1)
    resultSS->TD1 = 0;
  if (resultSS->PI_Net < 0)
    resultSS->PI_Net = 0;
  if (resultSS->PI_b < 0)
    resultSS->PI_b = 0;


  double *shannonIndexArray;

  shannonIndexArray = (double *) malloc (4 * sizeof (double));
  shannonIndex (list, n, &shannonIndexArray);
  resultSS->si1 = shannonIndexArray[0];
  resultSS->si2 = shannonIndexArray[1];
  resultSS->si3 = shannonIndexArray[2];
  resultSS->si4 = shannonIndexArray[3];

  free (shannonIndexArray);
  free(segwithin);
  free(freqdist);

  return resultSS;
}

/*
*/
void
PrintSumStatsArray (struct SumStat **SumStat_list, int numTaxonLocusPairs)
{
  int a;
  /* double MeanTAU, VarTAU, CV; */

  /* struct SumStat SumStat_list[numTaxonLocusPairs]; */
  /* WORK HERE, Mike I disabled sorting, this may need to be enabled later */
  //qsort (SumStat_list, numTaxonLocusPairs, sizeof (SumStat_list[0]), SS_comp);
  
  /****** NOTE ******
   *
   * (A) If new summary stat is added or the print order is
   *     changed, please modify the global: numStats and
   *     ssNameVect (top of this file).  numSumStats should be
   *     the number of summary statistics used for each taxon
   *     pair.
   *
   * (B) If new prior is added or the print order is changed,
   *     please modify numPriorColumns and priorNameVect.  For
   *     prior names, start with "PRI."
   *  
   * ORDER of names is important!
   */

  if (gPrintHeader)
    {
      int numPriorColumns = 0; /* Prior is not printed anymore */
      char priorNameVect[][MAX_LEN_COLUMN_NAME] =
	{ "PRI.Psi", "PRI.var.t", "PRI.E.t", "PRI.omega" };
      PrintHeader (priorNameVect, numPriorColumns, ssNameVect,
		   numSumStats, numTaxonLocusPairs);
      
      gPrintHeader = 0; /* lower the flag to print only once */
    }

  /* start to print sum stats */
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI_b);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI_w);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TW);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI_Net);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TD);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TDD);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI_w2);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->PI_w1);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TW2);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TW1);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TDD2);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->TDD1);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->si1);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->si2);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->si3);
  
  for (a = 0; a < numTaxonLocusPairs; a++)
    printf ("%lf\t", SumStat_list[a]->si4);
  
  printf ("\n");
    
  /*
  if ((fp=fopen("likeout1_21", "a+b")) ==NULL){
    fprintf(stderr,"Cannot open the file.\n");
    exit(1);
  }                                
  
  fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", numTauClassesHyper, VT, MeanTau, ObsvVARPi_Net, ObsvVARD, ObsvCV_pi_net_tW, ObsvEPi_Net );
  fclose (fp); 
  */
  
}

static void
PrintHeader (char priorNames[][MAX_LEN_COLUMN_NAME], int numPriors,
	     char sumStatNames[][MAX_LEN_COLUMN_NAME], int numSumStats,
	     int numTaxonPairs)
{
  int i, a;
  for (i = 0; i < numPriors; i++)
    {
      printf ("%s\t", priorNames[i]);
    }
  for (i = 0; i < numSumStats; i++)
    {
      for (a = 0; a < numTaxonPairs; a++)
	{
	  printf ("%s.%d\t", sumStatNames[i], a + 1);
	}
    }
  printf ("\n");
  return;
}

/*
 * Checks that sub population sample sizes n[] are reasonable.
 * Arguments:
 *   nsam:     number of total samples in the simulation
 *   npops:    number of subpopulations
 *   n[npops]: sub-population sample sizes
 *
 * Returns: 1 if  all subpop sample sizes are bet. 0 and nsam (ends exclusive)
 *          0 otherwise
 */
int
multiplepopssampledfrom (int nsam, int npops, int *n)	/*zzz i think this just tells the program to do Fst and that there is substructure zzz */
{
  int i, sum = 0;
  for (i = 0; i < npops; i++)
    {
      sum += n[i];

      if ((n[i] <= 0) || (n[i] >= nsam))
	return (0);
    }
  /* This function was checking only the first n[i], and returning 1
     if at least 1 element is 0 < n[i] < nsam.
     I don't think this is the intention, so I corrected it. Naoki
   */

  /* I have a feeling the additional check below may be good, too.  Naoki
     if (sum != nsam)
     return 0;
   */
  return (1);
}

int
pairwisediffc_w (int ss, int nsam, char **list, int np, int *n, int pop)
{
  int n1, n2, diffc = 0;
  int popi, startn = 0;
  int s;
  if (n[pop] > 1)
    {
      for (popi = 0; popi < pop; popi++)
	startn += n[popi];
      for (n1 = startn; n1 < (startn + n[pop] - 1); n1++)
	for (n2 = n1 + 1; n2 < (startn + n[pop]); n2++)
	  {
	    for (s = 0; s < ss; s++)
	      if (list[n1][s] != list[n2][s])
		diffc++;
	  }
    }
  /*printf("piW: %d\n", diffc);  test print */
  return (diffc);
}


int
pairwisediffc_b (int ss, int nsam, char **list, int np, int *n, int pop1,
		 int pop2)
{
  int n1, n2, diffc = 0;
  int popi, startn1, startn2;
  int s;
  if ((n[pop1] > 0) && (n[pop2] > 0))
    {
      startn1 = startn2 = 0;
      for (popi = 0; popi < pop1; popi++)
	startn1 += n[popi];
      for (popi = 0; popi < pop2; popi++)
	startn2 += n[popi];
      for (n1 = startn1; n1 < (startn1 + n[pop1]); n1++)
	for (n2 = startn2; n2 < (startn2 + n[pop2]); n2++)
	  {
	    for (s = 0; s < ss; s++)
	      if (list[n1][s] != list[n2][s])
		diffc++;
	    /*printf("diffc: %d", diffc);  test print */
	  }
    }
  /*printf("piB: %d\n", diffc);  test print */
  return (diffc);

}

/*yyy BELOW  Eli 05/15/06 yyy*/
/*
 * Arguments:
 *   segwReturnArr: results get returned here
 *   segsites: number of segSite in the entire data = # columns in list matrix
 *   list: sequence data array
 *   numSubpops: number of sub populations
 *   n: array with numSubpops elements.  i-th element is # seqs in i-th subpops
 *  
 * Assumes that memory is allocated to segwReturnArr (numSubpops * sizeof(int)).
 *
 * Side effect: segReturnArr get filled up by # seg sites within each subpops.
 */
void
segsubs (int *segwReturnArr, int segsites, char **list, int numSubpops, int *n)
{
  int i, count, npi, n0 = 0, ni, npops_gt1 = 0;

  /* initialize the result counter */
  memset(segwReturnArr, 0, numSubpops *  sizeof(int));

  for (npi = 0; npi < numSubpops; npi++)
    {
      if (n[npi] > 1)
	{
	  ++npops_gt1;
	  count = 0;
	  for (i = 0; i < segsites; i++)
	    {
	      for (ni = n0 + 1; ni < n0 + n[npi]; ni++)
		{
		  if (list[ni][i] != list[n0][i])
		    {
		      segwReturnArr[npi]++;
		      break;
		    }
		}
	    }
	}
      else
	{
	  segwReturnArr[npi] = 0;
	}
      n0 += n[npi];
    }
}

/*yyy ABOVE  Eli 05/15/06 yyy*/


/*
 * Calculate the pi (per gene and not per site) within sub populations.
 * Specify the index (0-offset) of sub population from which pi are calculated
 * in targetPop.
 *
 * If a negative value of targetPop is given, it calculate the pi within subpop
 * for each sub population and average is taken:  Let's say two subpops with
 * n0 and n1 samples, and pi's within each subpops are pi0 and pi1.
 * Weighted average is (C(n0,2) * n0 + C(n1,2) * n1) /(C(n0,2)+C(n1,2)),
 * where C(x,y) denote combinatorial: Choose y from x. 
 */
double
nucdiv_w (int nsam, int segsites, char **list, int np, int *n, int targetPop)
{
  int pop, pairwisediffc_w (int, int, char **, int, int *, int);
  int beginPop, endPop;
  double pi, nd;
  double num_comps;

  pi = 0.0;

  nd = nsam;
  num_comps = 0.;

  if (targetPop < 0)
    {
      beginPop = 0;
      endPop = np;
    }
  else
    {
      beginPop = targetPop;
      endPop = targetPop + 1;
    }

  for (pop = beginPop; pop < endPop; pop++)
    {
      pi += pairwisediffc_w (segsites, nsam, list, np, n, pop);
      num_comps += (double) n[pop] * ((double) n[pop] - 1.) / 2.;
    }

  pi /= num_comps;
  return (pi);
}


double
nucdiv_bw (int nsam, int segsites, char **list, int np, int *n)
{
  int pop1, pop2, pairwisediffc_b (int, int, char **, int, int *, int, int);
  double pi, nd;
  double num_comps;

  pi = 0.0;

  nd = nsam;
  num_comps = 0;
  for (pop1 = 0; pop1 < (np - 1); pop1++)
    for (pop2 = (pop1 + 1); pop2 < np; pop2++)
      {
	pi += pairwisediffc_b (segsites, nsam, list, np, n, pop1, pop2);
	/*printf("piB: %lf\n", pi);  test print */
	num_comps += (double) n[pop1] * (double) n[pop2];
      }
  pi /= num_comps;
  /*printf("piB-FINAL: %lf\n", pi);  test print */
  return (pi);
}


void
FrequencyDistrInfSites (int *freqdist, int n, int S, char **list)
{
  int i, oldfrequency (char, int, int, char **);
  for (i = 0; i < n; i++)
    freqdist[i] = 0;
  for (i = 0; i < S; i++)
    {
      freqdist[oldfrequency ('1', i, n, list)]++;	/* probably bogus for ACGT */
    }
}


void
FrequencyDistrQ (int *freqdist, int n, int S, char **list)
{
  int i, f, oldfrequency (char, int, int, char **);
  for (i = 0; i < n; i++)
    freqdist[i] = 0;
  for (i = 0; i < S; i++)
    {
      f = oldfrequency (list[0][i], i, n, list);
      freqdist[f < n / 2 + 0.0001 ? f : n - f]++;	/* probably bogus for ACGT */
    }
}


void
FuLi (D, F, n, S, list, pi)
     int n, S;
     char **list;
     double *D, *F, pi;
{
  int k, s, etae, oldfrequency (char, int, int, char **);
  double n1, S1, vD, uD, uF, vF, an, bn, cn;
  n1 = (double) n;
  S1 = (double) S;
  for (s = etae = 0; s < S; s++)
    if (oldfrequency ('1', s, n, list) == 1)
      etae++;
  for (k = 1, an = bn = 0.; k < n; k++)
    {
      an += 1. / (double) k;
      bn += 1. / (double) (k * k);
    }
  if (n == 2)
    cn = 1.;
  else
    cn = 2. * (n1 * an - 2. * (n1 - 1)) / ((n1 - 1) * (n1 - 2));
/*   printf("an:\t%9.7f\tbn:\t%9.7f\tcn:\t%9.7f\t",an,bn,cn); */
  vD = 1. + an * an / (bn + an * an) * (cn - (n1 + 1) / (n1 - 1));
  uD = an - 1. - vD;
  *D = S1 - an * etae;
  *D /= sqrt (uD * S1 + vD * S1 * S1);
  vF = cn + 2. * (n1 * n1 + n1 + 3) / (9. * n1 * (n1 - 1)) - 2. / (n1 - 1);
  vF /= an * an + bn;
  uF = 1. + (n1 + 1) / (3. * (n1 - 1)) - 4. * (n1 + 1) / 
    ((n1 - 1) * (n1 - 1)) * (an + 1. / n1 - 2. * n1 / (n1 + 1));
  uF /= an;
  uF -= vF;
/*   printf("vF:\t%9.7f\tuF:\t%9.7f\t",vF,uF); */
  *F = pi - etae;
  *F /= sqrt (uF * S1 + vF * S1 * S1);
}

double
tajddenominator (int n, int S, double p)
{
  int i;
  double n1, S1, a1, a2, b1, b2, e1, e2, denom;
  n1 = (double) n;
  S1 = (double) S;
  a1 = a2 = 0.;
  for (i = 1; i < n; i++)
    {
      a1 += 1. / i;
      a2 += 1. / (i * i);
    }
  b1 = (n1 + 1) / (3. * (n1 - 1));
  b2 = 2. * ((n1 * n1) + n1 + 3) / (9. * n1 * (n1 - 1));
  e1 = (b1 - 1. / a1) / a1;
  e2 = (b2 - (n1 + 2.) / (a1 * n1) + a2 / (a1 * a1)) / (a1 * a1 + a2);
  denom = sqrt (e1 * S1 + e2 * S1 * (S1 - 1));
  return (denom);
}

double
tajddenominator2 (int n, int S, double p)
{
  int i;
  double n1, S1, a1, a2, b1, b2, e1, e2, denom;
  n1 = (double) n;
  S1 = (double) S;
  a1 = a2 = 0.;
  for (i = 1; i < n; i++)
    {
      a1 += 1. / i;
      a2 += 1. / (i * i);
    }
  b1 = (n1 + 1) / (3. * (n1 - 1));
  b2 = 2. * ((n1 * n1) + n1 + 3) / (9. * n1 * (n1 - 1));
  e1 = (b1 - 1. / a1) / a1;
  e2 = (b2 - (n1 + 2.) / (a1 * n1) + a2 / (a1 * a1)) / (a1 * a1 + a2);
  denom = sqrt (e1 * S1 + e2 * S1 * (S1 - 1));
  return (denom);
}


double
thetaW (int n, int S)
{
  int i;
  double n1, S1, a1, theta;
  n1 = (double) n;
  S1 = (double) S;
  a1 = 0.;
  for (i = 1; i < n; i++)
    a1 += 1. / i;
  theta = S1 / a1;
  return (theta);
}


/* This calculate overall average pairwise differences (pi per gene)
 * ignoring sub population designation  */
double
nucdiv (int nsam, int segsites, char **list)
{
  int s;
  double pi = 0.0, denom;
  char dummy;

  for (s = 0; s < segsites; s++)
    {
      pi += frequency (dummy, s, nsam, list);
      /* frequency() returns the number of pair wise differences at site s 
       * from all pairwise comparison */
    }

  /* denom is  # of ways to choose a pair: nsam choose 2 */
  denom = nsam * (nsam - 1) / 2;
  return (pi / denom);
}


int
oldfrequency (char allele, int site, int nsam, char **list)
{
  int i, count = 0;
  for (i = 0; i < nsam; i++)
    count += (list[i][site] == allele ? 1 : 0);
  return (count);
}



/* 
 * Count the number of pairwise differences at the site.
 * nsam * (nsam - 1) / 2 pairs are compared.
 *
 * Arguments:
 *   base: ignored
 *   site: i-th segregating sites
 *   nsam: total number of samples
 *   list: character matrix of segregating sites
 *
 * Returns:  the number of pairwise differences at the site
 */
int
frequency (char base, int site, int nsam, char **list)
{
  char allele1;			/*7/27/04; Hickerson */
  int i, n, denom, count = 0;

  denom = 0;

  for (n = 0; n < nsam; n++)
    {
      allele1 = list[n][site];

      for (i = n; i < nsam; i++)
	{
	  if (list[i][site] != allele1)
	    count = count + 1;
	}

    }
  return (count);
}

int
frequencySING (char base, int site, int nsam, char **list)	/* in progress Hickerson 7/29/04 */
{
  char allele1;			/*7/27/04; Hickerson */
  int i, n, denom, singleton, count = 0;

  denom = 0;
  singleton = 0;
  for (n = 0; n < nsam; n++)
    {				/*7/27/04; Hickerson */
      allele1 = list[n][site];	/*7/27/04; Hickerson */
      count = 0;


      for (i = n; i < nsam; i++)
	{			/*7/27/04; Hickerson */
	  if (list[i][site] == allele1)
	    count = count;	/*7/27/04; Hickerson */
	  else
	    count = count + 1;	/*7/27/04; Hickerson */

	}
      if ((count = nsam - 1))
	{
	  singleton++;
	}
    }

  return (singleton);
}

/*  thetah - pi   */
/* 	double */
/* hfay( int nsam, int segsites, char **list) */
/* { */
/* 	int s, frequency( char, int, int, char**); */
/* 	double pi, p1, nd, nnm1  ; */

/* 	pi = 0.0 ; */

/* 	nd = nsam; */
/* 	nnm1 = nd/(nd-1.0) ; */
/*    	for( s = 0; s <segsites; s++){ */
/* 		p1 = frequency('1', s,nsam,list)/nd ; */
/* 		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ; */
/* 		} */
/* 	return( pi ) ; */
/* } */

/* Fay's theta_H  */
double
thetah (int nsam, int segsites, char **list)
{
  int s, oldfrequency (char, int, int, char **);
  double pi, p1, nd, nnm1;

  pi = 0.0;

  nd = nsam;
  nnm1 = nd / (nd - 1.0);
  for (s = 0; s < segsites; s++)
    {
      p1 = oldfrequency ('1', s, nsam, list);

      pi += p1 * p1;
    }
  return (pi * 2.0 / (nd * (nd - 1.0)));
}




int
segsub (int nsub, int segsites, char **list)
{
  int i, count = 0, c1;
  int oldfrequency (char, int, int, char **);

  for (i = 0; i < segsites; i++)
    {
      c1 = oldfrequency ('1', i, nsub, list);
      if ((c1 > 0) && (c1 < nsub))
	count++;
    }
  return (count);
}

void
shannonIndex (char **list, int *config, double **shannonIndexArray)
{
  int i, sizeOfSp1, sizeOfSp2, sizeAll, unit=1;
  double sHa1 = 0, sHa2 = 0, sHu = 0, sHua = 0, temp;
  hashtab_iter_t iHash;
  hashtab_t *subPop1, *subPop2, *pool;

  sizeOfSp1 = config[0];
  sizeOfSp2 = config[1];

  subPop1 = ht_init (sizeOfSp1, NULL);
  subPop2 = ht_init (sizeOfSp2, NULL);
  pool = ht_init (sizeOfSp1 + sizeOfSp2, NULL);

  if (subPop1 == NULL || subPop2 == NULL || pool == NULL) {
    fprintf(stderr, "ERROR: no memory in shannonIndex\n");
    exit(EXIT_FAILURE);
  }

  sizeOfSp1 = config[0];
  sizeOfSp2 = config[1];
  sizeAll = sizeOfSp1 + sizeOfSp2;

  /* insert allele-count pair into hashtables for subPop1
     key: allele as string, value: number of allele as int) */
  for (i = 0; i < sizeOfSp1; i++)
    {
      int *thisCnt = (int *) ht_search (subPop1, list[i], charCount (list[i]));
      if ( thisCnt == NULL) { /* new key */
	ht_insert (subPop1, list[i], charCount (list[i]), &unit, sizeof (int));
	ht_insert(pool,list[i], charCount (list[i]), &unit, sizeof (int));
      } else {/* this key have seen before */
	(*thisCnt)++;
	(*((int *) ht_search (pool, list[i], charCount (list[i])))) ++;
      }
    }

  /* For the subPop2 */
  for (i = 0; i < sizeOfSp2; i++)
    {
      int *thisCnt = (int*) ht_search (subPop2, list[sizeOfSp1 + i],
				       charCount (list[sizeOfSp1 + i]));
      
      if (thisCnt == NULL)
	ht_insert (subPop2, list[sizeOfSp1 + i], charCount(list[sizeOfSp1 + i]), 
		   &unit, sizeof (int));
      else
	(*thisCnt)++;

      /* deal with the pooled hashTbl */
      thisCnt = (int*) ht_search (pool, list[sizeOfSp1 + i],
				  charCount (list[sizeOfSp1 + i]));
      
      if (thisCnt == NULL)
	ht_insert (pool, list[sizeOfSp1 + i], charCount(list[sizeOfSp1 + i]), 
		   &unit, sizeof (int));
      else
	(*thisCnt)++;
      
    }

  // initialize hash table iterator, and calculate sHua1 
  for (ht_iter_init (subPop1, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / (double) sizeOfSp1;
      sHa1 += (-(log (temp) / log ((double) 2) * temp));
    }

  // calculate sHua2
  for (ht_iter_init (subPop2, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / sizeOfSp2;
      sHa2 += (-(log (temp) / log ((double) 2) * temp));
    }

  // calculate sHu
  for (ht_iter_init (pool, &iHash); iHash.key != NULL; ht_iter_inc (&iHash))
    {
      temp = (double) (*((int *) (iHash.value)));
      temp = temp / sizeAll;
      sHu += (-((log (temp) / log ((double) 2)) * temp));
    }

  // calculate sHua
  sHua = sHu - (sizeOfSp1 * sHa1 - sizeOfSp2 * sHa2)/sizeAll;

  // throw values into double array shannonIndexArray
  *(*shannonIndexArray + 0) = sHu, *(*shannonIndexArray + 1) =
    sHua, *(*shannonIndexArray + 2) = sHa1, *(*shannonIndexArray + 3) = sHa2;

  ht_destroy (pool);
  ht_destroy (subPop1);
  ht_destroy (subPop2);
}

/* 
 * Count the number of characters in a cstring (size of the char*)
 *
 * Argument:
 *   arr: the cstring whose number of characters to be counted
 *
 * Returns: the size of the string
 *
 */
int
charCount (char *arr)
{
  int k = 0;
  while (arr[k] != '\0')
    ++k;
  return k;
}				//int charCount(char*)
