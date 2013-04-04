#ifndef MS_PRIOR_H
#define MS_PRIOR_H

/* Default values used for interactive setup */
#define DEFAULT_THETA_SHAPE 1.0
#define DEFAULT_THETA_SCALE 0.001
#define DEFAULT_TAU_SHAPE  1.0
#define DEFAULT_TAU_SCALE  2.0
#define DEFAULT_BOTTLE_A 5.0
#define DEFAULT_BOTTLE_B 1.0
#define DEFAULT_MIGRATION_SHAPE -1.0
#define DEFAULT_MIGRATION_SCALE -1.0
#define DEFAULT_UPPER_REC  0.0
#define DEFAULT_ANCESTRAL_THETA_MULTIPLIER 1.0
#define DEFAULT_REPS  1000000
#define DEFAULT_MUT_FILE "SampleSize_MuModel_Vector"
#define DEFAULT_SUBPARAMCONSTRAIN "000000000"

#define MAX_FILENAME_LEN 1024
#define MAX_NAME_CHAR_LEN 1024
#define NUMBER_OF_CONPARAM 9

#include "hashtab.h"

// JRO - modified - 11/17/2011
typedef struct
{
  double thetaScale;
  double thetaShape;
  double tauScale;
  double tauShape;
  double bottleProportionShapeA;
  double bottleProportionShapeB;
  double migrationShape;
  double migrationScale;
  double upperRec;		/* upper limit of prior dist'n for recombination rate */
  double ancestralThetaMultiplier;
  unsigned long long reps;
  unsigned int numTaxonLocusPairs; /* total number of taxon:locus pairs */
  unsigned int numTaxonPairs;     /* number of unique taxon pairs */ 
  unsigned int numLoci;         /* number of unique loci */
  unsigned int numTauClasses;
  long prngSeed;
  char configFile[MAX_FILENAME_LEN];
  char scratchFile[MAX_FILENAME_LEN];
  unsigned int constrain;
  int printConf;                /* 1: -i opt print the config and exit (default: 0) */
  char subParamConstrain[NUMBER_OF_CONPARAM];
} runParameters;

extern runParameters gParam;

typedef struct
{
  char taxonName[MAX_NAME_CHAR_LEN];
  char locusName[MAX_NAME_CHAR_LEN];
  unsigned int taxonID;  /* integer (0-) representation of taxon pair name */
  unsigned int locusID;  /* integer (0-) representation of locus name */
  /* double thetaScaler; */
  double NScaler;     /* scaler of pop size, e.g. 0.25 for mtDNA of dioecy */
  double mutScaler;   /* mutation rate multiplier, relative to nuclear */
  unsigned int numPerTaxa;
  unsigned int sample[2];
  double tstv[2];
  double gamma;
  unsigned int seqLen;
  double freqA;
  double freqC;
  double freqG;
  double freqT;
  char filename[MAX_FILENAME_LEN];
} mutParameter;

typedef struct {
  int **tbl;
  unsigned int numTaxon;  /* number of rows (not 0-offset index) */
  unsigned int numLoci;   /* number of columns */
} lociTbl;

/* keeps track of sample size mutation model settings */
typedef struct
{
  mutParameter *data;
  int numElements;
  int numAllocated;
  int numSpecies;
  int numLoci;
  hashtab_t *taxonIDTbl;  /* taxon name string to -> unique ID hash table */
  hashtab_t *locusIDTbl;  /* locus name string to -> unique ID hash table */
  lociTbl *locTbl;
} mutParameterArray;

extern mutParameterArray gMutParam;

typedef struct
{
  double conTau;
  double conBottPop1;
  double conBottPop2;
  double conBottleTime;
  double conMig;		// migration rate
  double conTheta;
  double conN1;			//current population size1
  double conNanc;		// ancestral pupulation size
  double conRec;		// recombination rate
} constrainedParameter;

typedef struct
{
  constrainedParameter *conData;
  int conNumElements;
  int conNumAllocated;
} constrainedParameterArray;

extern constrainedParameterArray gConParam;

#endif
