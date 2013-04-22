#ifndef INTEGER_PARTITION_H
#define INTEGER_PARTITION_H

#include <gsl/gsl_rng.h>

int integerPartitionCumulativeSpectrum(int n, int arr[]);
int integerPartitionSpectrum(int n, int arr[]);
double integerPartitionCumulativeProbs(int n, double arr[]);
double integerPartitionProbs(int n, double arr[]);
int integerPartition(int n);
int drawIntegerPartitionCategory(gsl_rng * rng, int n);
int dirichletProcessDraw(gsl_rng * rng, int n, double alpha, int elements[]);
void generateIntegerPartitions(int n, int ip, int partitions[ip][n]);

#endif

