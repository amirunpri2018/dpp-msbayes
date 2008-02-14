/*
 * sumStatsVector.c
 *
 * Copyright (C) 2006  Michael Hickerson, Naoki Takebayashi
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

/*
  Change Log
  * Fri May 12 2006 Naoki Takebayashi <ffnt@uaf.edu>
  - Added ParseCommandLine().
  - tauequalizer of msprior and this program has to be synchronized.  So
    now this program takes upper bound of theta as an option (-T).
    This value enables us to calculate the tauequalizer.
  - msbayes.pl gets this value correctly from msprior, and pass it with -T.
  - If no -T is given, the default value of upperTheta which is hardcoded in
    msprior.h is used.
  - nadv used to be a commandline argument.  But I moved it to be an option.
    To pass nadv, we use -a N option now.
  - It will print the usage message with -h option.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <getopt.h> /* for getopt_long() */
#include <ctype.h> /* for isspace() */
#include "msprior.h"  /* MAX_FILENAME_LEN */

#define MAX_LNSZ 1000

/* posit[] was originally double, but I thougt int makes more sense, Naoki 
 * Change the typedef to double if I'm wrong */
typedef unsigned int tPositionOfSegSites;


runParameters gParam;
int gNadv = 0;
int gPrintHeader = 0;

/* function prototypes of external functions */
void printstats (int n, int S, char **, int subn, int npops, int *config,
		 double THETA, int isNumSegSitesConst, int Qbool, int Fst_bool,
		 double TAU, int count, int TAXAcount, int BasePairs,
		 int TauHyp, int NumTaxa);
int multiplepopssampledfrom (int nsam, int npops, int *config);
double thetaW (int, int), mu, N;

/* function prototypes for the functions in this file */
static char **cmatrix (int nsam, int len);
static int biggerlist (int nsam, unsigned nmax, char **list);
static void freeCMatrix(int nsam, char **list);
static int FindNumPopsAndSubpopSampleSizes (const char line[],
					    int **subPopSampleSize);
static void ReadInPositionOfSegSites(const char *line, 
				     tPositionOfSegSites *positionArray, 
				     int numSegSites);
char *FindFirstSpace(char *str);
char *RmLeadingSpaces(char *str);

static void ParseCommandLine(int argc, char *argv[]);
static int SetScratchFile(char *fName);

int main (int argc, char *argv[])
{
  int maxsites = 1000; /* max number of seg sites, used for size of data mat */
  int nsam, i, howmany, npops, *config;
  char **list, line[MAX_LNSZ], longline[32000], *mutscanline;
  char dum[100];
  FILE *pfin;
  tPositionOfSegSites *posit;
  double  theta, THETA, TAU, Tau1, Tau2;
  int segsites, count, TAXAcount,  BasePairs, NumTaxa;
  int TauHyp;
  int Fst_bool = 0, Qbool = 0;
  int isNumSegSitesConst = 0; /* 1 with -s, the number of segregating sites 
			       *    will be constant in each sample 
			       * 0 with -t, varies between samples
			       */
  char dum1[100], dum2[100], dum3[100], dum4[100], dum5[100];

  /* check the command line argument */
  /* if -T is given, the value from the option will override the default */
  gParam.upperTheta = DEFAULT_UPPER_THETA;  
  strncpy(gParam.scratchFile, "PARarray-E", MAX_FILENAME_LEN); /* set default scratch file */
  /* set gParam.upperTheta, gNadv (default 0) gParam.scratchFile */
  ParseCommandLine(argc, argv);  
  
  pfin = stdin;  /* read the data from STDIN */

  /* read in first line of output (command line options of msDQH) */
  fgets (line, MAX_LNSZ, pfin);

  /*
   * Get the following variables from the command line options
   *
   * nsam:      number of total samples
   * howmany:   how many simulations were run
   * THETA:     4 Ne mu used for the simulation
   * Tau1:      time of bottleneck event (going back in time)
   * Tau2:      time intvl between the bottleneck and separation events
   * TauHyp:    ???
   * BasePairs: sequence length
   * TAXAcount: sequential ID for each taxon pair (1 to # of taxon pairs)
   */  
  sscanf (line,
	  " %s %s %d %d %s %lf %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %lf %s %s %s %s %s %s %s %s %s %s %s %lf %s %s %s %u %s %s %s %u %s %s %s %d %s %s %s %u ",
	  dum, dum, &nsam, &howmany, dum, &THETA, dum, dum, dum, dum, dum, dum,
	  dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum, dum,
	  dum, dum, dum, dum, &Tau1, dum, dum, dum, dum, dum, dum, dum, dum,
	  dum, dum, dum, &Tau2, dum, dum5, dum, &TauHyp, dum4, dum, dum1,
	  &BasePairs, dum2, dum3, dum, &TAXAcount, dum, dum, dum, &NumTaxa);

  /* 
   * of course, I have to put a prior generator in the actual sample
   * generator for theta and tau down below for each count 
   */
  
  TAU = Tau1 + Tau2;           /* time of separation */
  TAU = TAU * (THETA * 2 / gParam.upperTheta);

  /* Find the theta or number of segregating sites from -t or -s */
  mutscanline = strstr (line, "-s");
  if (mutscanline != NULL)
    {
      /* number of segregating sites is constant */
      sscanf (mutscanline, " %d", &segsites);
      isNumSegSitesConst = 1;
      theta = thetaW (nsam, segsites);
    }
  else
    {
      mutscanline = strstr (line, "-t");
      if (mutscanline != NULL)
	sscanf (mutscanline, " %s", dum);
      else
	{
	  fprintf (stderr, "\nmutscanline problem -s or -t not found \n");
	  exit (1);
	}
      theta = atof (dum);

      /* -Q will tell transition transversion rate ratio and base freqs */
      if ((mutscanline = strstr (mutscanline, "-Q")) != NULL)
	Qbool = 1;
    }

  /* 
   * config become an array with npops elements, 
   * it contains subpop sample sizes
   */
  npops = FindNumPopsAndSubpopSampleSizes (line, &config);
  
  /* Checking if 0 < config[i] < nsam for all i */
  if ((npops > 1) && (multiplepopssampledfrom (nsam, npops, config)))
    Fst_bool = 1;

  /* prepare the storage for segregating sites data */
  if (isNumSegSitesConst)
    maxsites = segsites;

  list = cmatrix (nsam, maxsites + 1);
  posit = (tPositionOfSegSites *) calloc(maxsites,sizeof(tPositionOfSegSites));
  
  if (list == NULL || posit == NULL) {
    fprintf(stderr, "No mem for segregating sites data\n");
    exit(EXIT_FAILURE);
  }

  /* Start to process the data */
  count = 0;
  while (howmany - count++)
    {
      /* The line after "//" is the beginning of simulation data */
      while (strcmp(line, "//\n") != 0)
	fgets (line, MAX_LNSZ, pfin);
      
      /* Number of segregating sites line */
      fgets (line, MAX_LNSZ, pfin);
      if (! isNumSegSitesConst) {
	sscanf (line, "segsites: %d\n", &segsites);
	
	if (segsites >= maxsites)      /* readjust the size of data matrix */
	  {
	    maxsites = segsites + 10;  /* extra 10 elements */
	    posit = (tPositionOfSegSites *) 
	      realloc(posit, maxsites * sizeof (tPositionOfSegSites));
	    /*printf("PRE %d %d %d\n", segsites, maxsites, nsam);*/ 
	    if (posit == NULL || biggerlist(nsam, maxsites, list) != 0 ) {
	      fprintf(stderr,
		      "Not enough memory for reallocating char matrix\n");
	      exit(EXIT_FAILURE);
	    }
	  }
      }

      /* get rid of base frequency line */
      if (Qbool)
	{
	  fgets (line, MAX_LNSZ, pfin);
	  sscanf (line, "freqACGT: %s %s %s %s", dum, dum, dum, dum);
	}
      
      if (segsites > 0)
	{
	  /* read in position of segregating sites */
	  fgets (longline, 32000, pfin);

	  /* posit array initialized */
	  ReadInPositionOfSegSites(longline, posit, segsites);

	  /* list[][] get initialized with character states */
	  for (i = 0; i < nsam; i++)
	    fscanf (pfin, " %s", list[i]);

	}
      /* what do we do if segsites = 0?, Naoki */

      /* analyse sample ( do stuff with segsites and list) */
      printstats (nsam, segsites, list, gNadv, npops, config, THETA, 
		  isNumSegSitesConst, Qbool, Fst_bool, TAU, count, TAXAcount, 
		  BasePairs, TauHyp, NumTaxa);

      /* perhaps put acceptance/rejection thing here */
      
    }

  free (config);
  free(posit);
  freeCMatrix(nsam, list);
  return (0);
}


/*
 * Take a character string line[] as the argument.  Find -m (migration
 * rate) or -D (demographic history?) from this line[], and extract
 * the information about sub-population sample sizes.
 *
 * The line[] is something like this:
 *
./msDQH 15 1 -t 22.052669 -Q 4.765400 0.207500 0.231800 0.215000 0.345700 -H 999.000000 -r 137.043900 639 -D 5 2 2 13 0 I 0.000000 1.723360 1.723360 0.276640 0.276640 2.522721 2 1 0 0 1 0 I 0.000000 Nc 0.352105 0.673049 4.170096 1 Nc 0.215405 1.1 1 Nc 0.215405 639 1 Nc 0.215405 2
 * 
 *
 * Returned value: numSubPops
 *     The number of subpopulations (1st number after these options) 
 * 
 * Side effect:
 *   Mem is allocated to the pointer *subpopSampleSize (numSubPops elements),
 *   and this array contains the sample sizes of sub-populations.
 *
 */
static int 
FindNumPopsAndSubpopSampleSizes (const char line[], int **subPopSampleSize){
  char * mscanline, *Dscanline, *charPtr;
  int numSubPops, i;
  int *arrayPtr;

  mscanline = strstr (line, "-m");
  Dscanline = strstr (line, "-D");

  if ((mscanline == NULL) && (Dscanline == NULL)) {
    *subPopSampleSize = (int *) malloc(sizeof(int));
    if (*subPopSampleSize == NULL) {
      fprintf(stderr, "No memory in FindNumPopsAndSubpopSampleSizes()\n");
      exit(EXIT_FAILURE);
    }
    return 1;  /* no -D nor -m, so only 1 population */
  }

  /* set the char pointer to the beginning of arguments for the option */
  if (mscanline != NULL ) {
    if (Dscanline != NULL) /* Both -m & -D specified, can be a problem, Naoki*/
      fprintf(stderr, "WARN: both -m and -D are specified, ignoring -D\n");
    
    charPtr = mscanline;

  } else if (Dscanline != NULL) {

    charPtr = Dscanline;
  }

  /* the 2nd argument after -D contains the number of sub populations. */
  for (i = 0; i <2; i++) {
    charPtr = FindFirstSpace(charPtr);
    if(! charPtr) {
      fprintf(stderr, "ERROR: no arguments are given after -m or -D");
      exit(EXIT_FAILURE);
    }
    charPtr = RmLeadingSpaces(charPtr);
  }
  numSubPops = (int) strtol (charPtr, &charPtr, 10);
  
  /* subPopSampleSize is an array storing the sub-population sample sizes */
  arrayPtr = (int *) calloc (numSubPops, sizeof (int));
  /* This sucks up the subpop sample sizes into arrayPtr[] */
  for (i = 0; i < numSubPops; i++) {
    arrayPtr[i] = (int) strtol (charPtr, &charPtr, 10);
    
    if(errno == ERANGE) {
      fprintf(stderr, "WARN: out of range in population sizes. "
	      "Continuing... But the results are probably wrong\n");
    }
    if (errno == EINVAL) {
      fprintf(stderr, "ERROR: The arguments to -m options is weird:\n"
	      "%s\n", line);
      exit(EXIT_FAILURE);
    }
  }
  
  *subPopSampleSize = arrayPtr;  /* returning the address of the array */
  return numSubPops;
}


/*
 * Process a line containing positions of segregating sites, and
 * assign the values to an array
 *
 * Arguments:
 *   line: A character string which contains the position of segregating sites
 *         This line starts with "positions:" and contains integers 
 *         delimited by spaces.  There should be numSegSites integers.
 *
 *   positionArray: The values extracted from line will be assigned to this
 *                  array.  Assumes correct size mem is allocated
 *   numSegSites: number of segregating sites.  The size of positionArray
 *                should be equal to (or greater than) this value.
 *
 * Returns: nothing
 *
 * Side effects: values will be assigned to positionArray
 *     
 */
static void
ReadInPositionOfSegSites(const char *line, tPositionOfSegSites *positionArray, 
			 int numSegSites)
{
  int i;
  char *charPtr, *tempPtr;
  
  charPtr = strstr(line, "positions:");
  if (charPtr == NULL) {
    fprintf(stderr, "ERROR: positions line not found, ignoring\n");
    exit(EXIT_FAILURE);
  }
  
  charPtr = FindFirstSpace(charPtr);
  if (charPtr == NULL) {
    fprintf(stderr, "ERROR: positions line doesn't contain a space\n");
    exit(EXIT_FAILURE);
  }
  
  for (i = 0; i < numSegSites; i++) {
    positionArray[i] = (int) strtol (charPtr, &tempPtr, 10);
    if (charPtr == tempPtr) {
      fprintf(stderr, "ERROR: reached to the end, while processing "
	      "positions line\n");
      exit(EXIT_FAILURE);
    }
    charPtr = tempPtr;
  }
  return;
}

/*****************  Character matrix functions **********************/
/* allocates space for gametes (character strings) */
static char ** cmatrix (int nsam, int len)
{
  int i;
  char **m;
  if (!(m = (char **) malloc ((unsigned) (nsam * sizeof (char *))))) {
    perror ("alloc error in cmatrix");
    return NULL;
  }
  for (i = 0; i < nsam; i++)
    {
      if (!(m[i] = (char *) malloc ((unsigned) (len * sizeof (char))))) {
	perror ("alloc error in cmatric. 2");
	return NULL;
      }
    }
  return (m);
}

/* 
 * Arguments:
 *   list[][]: a matrix storing segregating sites data, allocated by cmatrix()
 *   nsam:     number of rows in list[][]
 *   nmax:     number of columns
 * 
 * Side effect:
 *   grow the size (character string length) to nmax.
 *
 * Returns 0 on success
 *        -1 on failure
 */   
static int biggerlist (int nsam, unsigned nmax, char **list)
{
  int i;

  for (i = 0; i < nsam; i++)
    {
      list[i] = (char *) realloc (list[i], nmax * sizeof (char));
      if (list[i] == NULL) {
	perror ("realloc error. in biggerlist()");
	return -1;
      }
    }
  return 0;
}

static void freeCMatrix(int nsam, char **list) {
  int i;
  
  for (i = 0; i < nsam; i++)
    free(list[i]);

  free(list);
  return;
}

/* Find the first space character and return the pointer to it 
 * Returns NULL if no space is found
 */
char *FindFirstSpace(char *str) {
  char *cPtr;
  if (str) {
    for (cPtr = str; *cPtr && ! (isspace(*cPtr)); cPtr++) {
      ;
    }
    if (cPtr != str + strlen(str)) {
      return cPtr;
    }
  }
  
  return NULL;
}

/* Macro to move string from s to d */
#define strMove(d,s) memmove(d,s,strlen(s)+1)

/*
 *  Remove leading white spaces from a string
 */
char *RmLeadingSpaces(char *str)
{
  char *obuf;

  if (str)
    {
      for (obuf = str; *obuf && isspace(*obuf); ++obuf)
        ;
      if (str != obuf)
        strMove(str, obuf);
    }
  return str;
}

/*
 * Print out the usage, and exit
 */
static void PrintUsage(char *progname)
{
  char *p;

  /*
   * Strip the path from the program name for when we
   * print the usage.
   */
  p = strrchr(progname, '/');
  if(p)
    p++;
  else
    p = progname;

  fprintf(stderr,
          "\nUsage: %s [--help] [--header] [--upper_theta N] [--nadv N] [--tempFile fileName] [< output_line_of_msDQH]\n\n"
          "        help: Print this usage function (-h)\n"
	  "      header: Print column header (-H)\n"
          " upper_theta: Specify a upper limit of prior distribution for "
	                 "theta (-T)\n\n"
	  "        nadv: Specify nadv (-a)\n"
	  "     tmpFile: Specify a filename which can be used to store temporary data (-t)\n"
          "stdin is used to read in a single line of msDQH output "
	  "(output_line_of_msDQH)"
          "\n\n", p);
  exit(EXIT_FAILURE);
}


static struct option sim_opts[] = {
  { "help", 0, NULL, 'h'},  /* list options */
  { "header",0,NULL, 'H'},
  { "upper_theta", 1, NULL, 'T'},
  { "nadv", 1, NULL, 'a'},
  { "tempFile", 1, NULL, 't'},
  { NULL, 0, NULL, 0}
};

static void ParseCommandLine(int argc, char *argv[])
{
  while(1) {
    int opt, rc;
    int opt_index;
    
    opt = getopt_long(argc, argv, "hHT:a:t:", sim_opts, &opt_index);
    if(opt < 0)
      break;
    
    switch(opt) {
    case 'h':   /* Print usage and exit */
      PrintUsage(argv[0]); /* This function will exit */
      break;
    case 'H':   /* Print header */
      gPrintHeader = 1;
      break;
    case 'T':   /* Specify the upperTheta */
      if (!optarg) {
	fprintf(stderr, "Must select upper bound of prior distribution "
		"of theta\n");
	PrintUsage(argv[0]);
      }
      gParam.upperTheta = strtod(optarg, NULL);
      if(errno || (gParam.upperTheta < 0)) {
	fprintf(stderr, "Invalid upper theta: %s\n", optarg);
	PrintUsage(argv[0]);
      }
      break;
    case 'a':  /* specify nadv */
      if (!optarg) {
	fprintf(stderr, "Must select nadv value "
		"of theta\n");
	PrintUsage(argv[0]);
      }
      gNadv = strtol(optarg, NULL, 10);
      if(errno || (gNadv < 0)) {
	fprintf(stderr, "Invalid nadv: %s\n", optarg);
	PrintUsage(argv[0]);
      }
      break;
    case 't':
      if (!optarg) {
	fprintf(stderr, "Must specify filename used for temporary storage file\n");
	PrintUsage(argv[0]);
      }
      rc=SetScratchFile(optarg);
      if(rc<0) {
	fprintf(stderr, "Can't set the temporary scratch file as '%s'\n", optarg);
	PrintUsage(argv[0]);
      }
      break;
    default:
      PrintUsage(argv[0]); /* This function will exit */
      break;
    }
  }
}


static int SetScratchFile(char *fName) 
{
  int len;

  if(!fName || !fName[0])
    return -1;

  len = strlen(fName);
  if((len + 1) >= MAX_FILENAME_LEN)
    return -1;
  
  strncpy(gParam.scratchFile, fName, MAX_FILENAME_LEN);

  return 0;
}
