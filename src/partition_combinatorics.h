#ifndef PARTITION_COMBINATORICS_H
#define PARTITION_COMBINATORICS_H

#include <gsl/gsl_rng.h>

int integerPartitionCumulativeSpectrum(int n, int arr[]);
int integerPartitionSpectrum(int n, int arr[]);
double integerPartitionCumulativeProbs(int n, double arr[]);
double integerPartitionProbs(int n, double arr[]);
int integerPartition(int n);
int drawIntegerPartitionCategory(const gsl_rng * rng, int n);
int dirichletProcessDraw(const gsl_rng * rng, int n, double alpha, int elements[]);
void generateIntegerPartitions(int n, int ip, int **partitions);

#endif

