/**
 * \file rejection.c
 * \brief A collection of functions for Euclidean distance rejection.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct vect_ {
    double * v;
    int length;
} vect;

vect init_vect(int length) {
    vect v;
    v.length = length;
    if ((v.v = (typeof(*v.v) *) calloc(v.length, sizeof(*v.v))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    return v;
}

double get_euclidean_distance(vect v1, vect v2) {
    if (v1.length != v2.length) {
        fprintf(stderr, "ERROR: get_euclidean_distance: arrays must be of "
                "equal length\n");
        exit(1);
    }
    double sum_of_squared_diffs = 0.0;
    int i;
    for (i = 0; i < v1.length; i++) {
        sum_of_squared_diffs += pow((v1.v[i] - v2.v[i]), 2);
    }
    return sqrt(sum_of_squared_diffs);
}

void standardize_vector(vect v, vect means, vect std_devs) {
    if ((v.length != means.length) || (v.length != std_devs.length)) {
        fprintf(stderr, "ERROR: standardize_vector: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < v.length; i++) {
        v.v[i] = (v.v[i] - means.v[i]) / std_devs.v[i];
    }
}

int main(int argc, char **argv){
    int i, l = 4;
    vect v1, v2, means, std_devs;
    v1 = init_vect(l);
    v2 = init_vect(l);
    means = init_vect(l);
    std_devs = init_vect(l);
    for (i = 0; i < l; i++) {
        v1.v[i] = 2;
        v2.v[i] = 4;
        means.v[i] = 3;
        std_devs.v[i] = 2;
    }
    double d1 = get_euclidean_distance(v1, v2);
    standardize_vector(v1, means, std_devs);
    standardize_vector(v2, means, std_devs);
    double d2 = get_euclidean_distance(v1, v2);
    fprintf(stdout, "%f\n%f\n", d1, d2);
    vect v3 = init_vect(5);
    double d3 = get_euclidean_distance(v1, v3);
    return 0;
}

