#include "partitionCombinatorics.h"

int integerPartitionCumulativeSpectrum(int n, int arr[]) {
    int i, j;
    int table[n+1][n+1];
    // initialize table with base cases
    for (i = 0; i <= n; i++)
        table[i][0] = 0;
    for (i = 1; i <= n; i++)
        table[0][i] = 1;
    // populate table
    for (i = 1 ; i <= n; i++) {
        for (j = 1; j <= n; j++) {
            if (i - j < 0){
                table[i][j] = table[i][j-1];
                continue;
            }
            table[i][j] = table[i][j-1] + table[i-j][j];
        }
    }
    for (i = 0; i < n; i++) {
        arr[i] = table[n][i+1];
    }
    return table[n][n];
}

int integerPartitionSpectrum(int n, int arr[]) {
    int i;
    int a[n];
    int ip = integerPartitionCumulativeSpectrum(n, a);
    arr[0] = 1;
    for (i = 1; i < n; i++) {
        arr[i] = a[i] - a[i-1];
    }
    return ip;
}
    
double integerPartitionCumulativeProbs(int n, double arr[]) {
    int i;
    int a[n];
    int ip = integerPartitionCumulativeSpectrum(n, a);
    for (i = 0; i < n; i++) {
        arr[i] = a[i] / ((double) ip);
    }
    return arr[n-1];
}

double integerPartitionProbs(int n, double arr[]) {
    int i;
    int a[n];
    int ip = integerPartitionSpectrum(n, a);
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        arr[i] = a[i] / ((double) ip);
        sum += arr[i];
    }
    return sum;
}

double integerPartition(int n) {
    int arr[n];
    int ip = integerPartitionCumulativeSpectrum(n, arr);
    return ip;
}

int drawIntegerPartitionCategory(gsl_rng * rng, int n) {
    double cumulativeProbs[n];
    double totalProb = integerPartitionCumulativeProbs(n, cumulativeProbs);
    double r = gsl_rng_uniform(rng);
    int i;
    for (i = 0; i < n; i++) {
        if (r < cumulativeProbs[i]) {
            return (i+1);
        }
    }
    return -1;
}

int dirichletProcessDraw(gsl_rng * rng, int n, double alpha, int elements[]) {
    elements[0] = 0;
    int subSetCounts[n];
    subSetCounts[0] = 1;
    int numSubSets = 1;
    double subSetProb, newSubSetProb, u;
    int i, j;
    for (i = 1; i < n; i++) {
        newSubSetProb = (alpha / (alpha + (double)i));
        u = gsl_rng_uniform(rng);
        u -= newSubSetProb;
        if (u < 0.0) {
            elements[i] = numSubSets;
            subSetCounts[numSubSets] = 1;
            numSubSets += 1;
            continue;
        }
        for (j = 0; j < numSubSets; j++) {
            subSetProb = ((double)subSetCounts[j] / (alpha + (double)i));
            u -= subSetProb;
            if (u < 0.0) {
                elements[i] = j;
                subSetCounts[j] += 1;
                break;
            }
        }
        if (u > 0.0) {
            elements[i] = numSubSets-1;
            subSetCounts[numSubSets-1] += 1;
        }
    }
    return numSubSets;
}

