#ifndef MS_PRIOR_H
#define MS_PRIOR_H

/* Default values used for interactive setup */
#define DEFAULT_LOWER_THETA 0.5
#define DEFAULT_UPPER_THETA 40.5
#define DEFAULT_UPPER_TAU  10.0
#define DEFAULT_UPPER_MIG  0.0 /* when set 0, it is a constant */
#define DEFAULT_UPPER_REC  0.0
#define DEFAULT_UPPER_ANC_POPSIZE 0.5
#define DEFAULT_NUM_LOCI 1;
#define DEFAULT_REPS  500000
#define DEFAULT_MUT_FILE "SampleSize_MuModel_Vector"
#define DEFAULT_PRIOR_OUT_FILE "msBayesPriorOut.csv"
#define DEFAULT_SUBPARAMCONSTRAIN "000000000"

#define MAX_FILENAME_LEN 255
#define NUMBER_OF_CONPARAM 9

typedef struct {
  double upperTheta;/* upper limit of prior dist'n for theta */
  double lowerTheta;/* lower limit of prior dist'n for theta */
  double upperTau;  /* upper limit of prior dist'n for time of separation */
  double upperMig;  /* upper limit of prior dist'n for migration rate */
  double upperRec;  /* upper limit of prior dist'n for recombination rate */
  double upperAncPopSize; /* upper limit of prior dist'n for ancestral pop size */
  unsigned long long reps;
  unsigned int numTaxaPair;
  unsigned int numTauClasses;
  unsigned int numLoci;
  long prngSeed;
  char configFile[MAX_FILENAME_LEN];
  char priorOutFile[MAX_FILENAME_LEN];
  char scratchFile[MAX_FILENAME_LEN];
  unsigned int constrain;
  char subParamConstrain[NUMBER_OF_CONPARAM];
} runParameters;

extern runParameters gParam;

typedef struct {
  unsigned int numPerTaxa;
  unsigned int sample[2];
  double tstv[2];
  double gamma;
  unsigned int seqLen;
  double freqA;
  double freqC;
  double freqG;
  double freqT;
} mutParameter;

typedef struct{
  mutParameter *data;
  int numElements;
  int numAllocated;
} mutParameterArray;

extern mutParameterArray gMutParam;

typedef struct{
  double conTau;
  double conBottPop1;
  double conBottPop2;
  double conBottleTime;
  double conMig; // migration rate
  double conTheta; 
  double conN1; //current population size1
  double conNanc; // ancestral pupulation size
  double conRec; // recombination rate
} constrainedParameter;

typedef struct{
  constrainedParameter *conData;
  int conNumElements;
  int conNumAllocated;
} constrainedParameterArray;

extern constrainedParameterArray gConParam;

#endif
