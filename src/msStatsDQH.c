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
/* #include <float.h> */
/* #if defined (__SVR4) && defined (__sun) */
/*   int isinf(double x) { return !finite(x) && x==x; } */
/* #endif */

double nucdiv (int, int, char **);
double nucdiv_w (int, int, char **, int, int *);
double nucdiv_bw (int, int, char **, int, int *);
double tajddenominator (int, int, double);
double thetaW (int, int);
double thetah (int, int, char **);
void FuLi (double *D, double *F, int, int, char **, double pi);
void FrequencyDistrInfSites (int *freqdist, int nsam, int segsites,
			     char **list);
void FrequencyDistrQ (int *freqdist, int nsam, int segsites, char **list);
int segsub (int nsub, int segsites, char **list);
int multiplepopssampledfrom (int nsam, int npops, int *n);	/*zzz n[] is called config[] in samples.c zzz */
static int SS_comp (const void *, const void *);
#if 0
static int compare_doubles (const void *a, const void *b);
#endif

struct SumStat
{
  double PI_Net;
  double PI_b;
  double PI_w;
  double PI;
  double TD;
  double TW;
  double TDD;
};

#if 0  /* commented out since it is not used */
/*
 * used for qsort to compare two doubles.
 * Takes two pointers to double as the arguments.
 *
 * Returns: 1  if a > b
 *          0  if a == b
 *         -1  if a < b
 */

static int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  
  return (*da > *db) - (*da < *db);
}
#endif


static int SS_comp (const void *p1, const void *p2)
{
  const struct SumStat *sp1 = (struct SumStat *) p1;
  const struct SumStat *sp2 = (struct SumStat *) p2;

  return ((sp1->PI_Net) > (sp2->PI_Net)) - ((sp1->PI_Net) <
					    (sp2->PI_Net));
  /*return  ( sp1->PI_b) - ( sp2->PI_b); */
}


void
printstats (int nsam, int segsites, char **list, int nsub, int npops, int *n,
	    double THETA, int Sbool, int Qbool, int Fst_bool, double TAU,
	    int count, int TAXAcount, int BasePairs, int TauHyp, int NumTaxa)
{
  int  STATLOAD, SSLOAD;
  double tW, pi, D,  FuLiD, FuLiF, CV, TDen;
  double pi_w=-1, pi_b=-1, Fst, Nm, Pi_Net = -1;/*zz7z Hickerson 7/29/04 zzz */
  int *freqdist, nsegsub=-1, CC, a;
  double NSEGSUB[NumTaxa];
  double TW[NumTaxa], PI[NumTaxa], PI_w[NumTaxa], PI_b[NumTaxa], 
    PI_Net[NumTaxa], TD[NumTaxa], tau[NumTaxa], TDD[NumTaxa];
  double MeanTAU, VarTAU;
  /* double h, th, ObsvVARD, ObsvVARPi_Net, ObsvEPi_Net, ObsvCV_pi_net_tW; */

  FILE *fp;

  struct SumStat SumStat_list[NumTaxa];



  freqdist = (int *) malloc (nsam * sizeof (int));
  if (Qbool)
    {
      FrequencyDistrQ (freqdist, nsam, segsites, list);
    }
  else
    {
      FrequencyDistrInfSites (freqdist, nsam, segsites, list);
    }
  if (Sbool)
    tW = THETA;
  else
    tW = (double) thetaW (nsam, segsites);
  pi = (double) nucdiv (nsam, segsites, list);
  /*  printf("%d\n", BasePairs); */


  /*  printf("%d\t%d\n", npops, n); */
  if (Fst_bool)
    {
      pi_w = (double) nucdiv_w (nsam, segsites, list, npops, n) / BasePairs;
      pi_b = (double) nucdiv_bw (nsam, segsites, list, npops, n) / BasePairs;
      Fst = 1. - pi_w / pi_b;
      Pi_Net = (double) pi_b - pi_w;
      if (Fst < 0)
	{
	  Fst = 0.;
	  Nm = -1.;
	}
      else
	Nm = (1. / Fst - 1.) / 4.;
    }				/*zzz Hickerson 7/29/04 zzz */
/*   th = thetah(nsam, segsites, list) ; */
  D = (pi - tW) / tajddenominator (nsam, segsites, pi);
  TDen = tajddenominator (nsam, segsites, pi);
/*  FuLi(&FuLiD,&FuLiF,nsam,segsites,list,pi);*/
/*   h = pi-th ; */
  if (nsub > 0)
    nsegsub = segsub (nsub, segsites, list);


  tW = tW / BasePairs;
  pi = pi / BasePairs;

  if (segsites < 1)
    D = 0, Fst = 0, Pi_Net = 0, FuLiD = 0, FuLiF = 0, tW = 0, pi = 0, pi_b =
      0, pi_w = 0, TDen = 0;
  if (Pi_Net < 0)
    Pi_Net = 0;
  if (pi_b < 0)
    pi_b = 0;


  CC = TAXAcount - 1;
 
  
  if ((fp = fopen ("PARarray-E", "a+b")) == NULL)
    {
      fprintf (stderr, "Cannot open the file.\n");
      exit (1);
    }


  /*fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", pi_b, FuLiD, FuLiF, tW, pi, segsites, pi_w, Fst, Pi_Net, D, TAU, TDen);
     fclose (fp); */

  fprintf (fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", pi_b, tW, pi, pi_w,
	   Pi_Net, D, TDen, TAU);
  fclose (fp);


  if (nsub > 0)
    NSEGSUB[CC] = nsegsub;

  if (TAXAcount == NumTaxa)		
    {
      if ((fp = fopen ("PARarray-E", "rb")) == NULL)
	{
	  fprintf (stderr, "Cannot open the file.\n");
	  exit (1);
	}
      STATLOAD = NumTaxa;
      SSLOAD = 11;
      for (a = 0; a < STATLOAD; a++)
	fscanf (fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &PI_b[a],
		&TW[a], &PI[a], &PI_w[a], &PI_Net[a], &TD[a], &TDD[a],
		&tau[a]);

      fclose (fp);

      for (a = 0; a < STATLOAD; a++)
	{
	  SumStat_list[a].PI_Net = PI_Net[a];
	  SumStat_list[a].PI_b = PI_b[a];
	  SumStat_list[a].PI_w = PI_w[a];
	  SumStat_list[a].PI = PI[a];
	  SumStat_list[a].TD = TD[a];
	  SumStat_list[a].TW = TW[a];
	  SumStat_list[a].TDD = TDD[a];
	}
 
      qsort (SumStat_list, NumTaxa, sizeof SumStat_list[0], SS_comp);

      {
	VarTAU = gsl_stats_variance (tau, 1, NumTaxa);
	MeanTAU = gsl_stats_mean (tau, 1, NumTaxa);
	CV = VarTAU / MeanTAU;


	printf ("%d\t%lf\t%lf\t%lf\t", TauHyp, VarTAU, MeanTAU, CV);
	/*for (a=0;a<STATLOAD;a++)printf("%lf\t", SumStat_list[a].PI_Net);
	   printf("\n"); */

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].PI_b);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].PI_w);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].PI);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].TW);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].PI_Net);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].TD);

	for (a = 0; a < STATLOAD; a++)
	  printf ("%lf\t", SumStat_list[a].TDD);
	printf ("\n");


	/*    if ((fp=fopen("likeout1_21", "a+b")) ==NULL){
	   fprintf(stderr,"Cannot open the file.\n");
	   exit(1);
	   }                                

	   fprintf(fp, "%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\n", TauHyp, VT, MeanTau, ObsvVARPi_Net, ObsvVARD, ObsvCV_pi_net_tW, ObsvEPi_Net );
	   fclose (fp); */

      }

      remove ("PARarray-E");

      /*remove arrays */


    }
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
  for (i = 0; i < npops; i++) {
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

double
nucdiv_w (int nsam, int segsites, char **list, int np, int *n)
{
  int pop, pairwisediffc_w (int, int, char **, int, int *, int);
  double pi, nd;
  double num_comps;

  pi = 0.0;

  nd = nsam;
  num_comps = 0.;
  for (pop = 0; pop < np; pop++)
    {
      pi += pairwisediffc_w (segsites, nsam, list, np, n, pop);
      /* printf("piW: %lf\n", pi);  test print */
      num_comps += (double) n[pop] * ((double) n[pop] - 1.) / 2.;
    }
  pi /= num_comps;
  /*printf("piW-FINAL: %lf\n", pi);  test print */
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
  uF =
    1. + (n1 + 1) / (3. * (n1 - 1)) - 4. * (n1 +
					    1) / ((n1 - 1) * (n1 - 1)) * (an +
									  1. /
									  n1 -
									  2. *
									  n1 /
									  (n1
									   +
									   1));
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


/* Calcuate the average pairwise differences */
double
nucdiv (int nsam, int segsites, char **list)
{
  int s, frequency (char, int, int, char **);	/*7/27/04; Hickerson */
  double pi,  denom;
  char dummy = '?';
  pi = 0.0;
  for (s = 0; s < segsites; s++)
    {
      /* 
       * frequency() returns the number of pair wise differences
       * at site s from all pairwise comparison
       */
      pi += frequency (dummy, s, nsam, list);	/*7/27/04; Hickerson */
    }
  denom = nsam * (nsam - 1) / 2;   /* # of ways to choose a pair:
				      nsam choose 2 */
  pi = pi / denom;
  /*printf("piFINAL: %lf\n", pi);  test print */

  return (pi);
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
    {				/*7/27/04; Hickerson */
      allele1 = list[n][site];	/*7/27/04; Hickerson */

      /*printf("n: %d\t site; %d\t nsam: %d\t allele1: %c\n", n, site, nsam, allele1);  test print */


      for (i = n; i < nsam; i++)
	{			/*7/27/04; Hickerson */
	  if (list[i][site] != allele1)
	    count = count + 1;	/*7/27/04; Hickerson */
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
      if ( (count = nsam - 1) ) {
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
