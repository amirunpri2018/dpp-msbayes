/**
 * @file        partition_combinatorics.c
 * @authors     Jamie Oaks
 * @package     msBayes
 * @brief       A collection of functions for partition combinatorics.
 * @copyright   Copyright (C) 2013 Jamie Oaks.
 *   This file is part of msBayes.  msBayes is free software; you can
 *   redistribute it and/or modify it under the terms of the GNU General Public
 *   License as published by the Free Software Foundation; either version 2 of
 *   the License, or (at your option) any later version.
 * 
 *   msBayes is distributed in the hope that it will be useful, but WITHOUT ANY
 *   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *   FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *   details.
 * 
 *   You should have received a copy of the GNU General Public License along
 *   with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "partition_combinatorics.h"

int integer_partition_cumulative_spectrum(int n, i_array * spectrum_dest) {
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
    (*spectrum_dest).length = 0;
    for (i = 0; i < n; i++) {
        append_i_array(spectrum_dest, table[n][i+1]);
    }
    return table[n][n];
}

int integer_partition_spectrum(int n, i_array * spectrum_dest) {
    int i;
    i_array v;
    v = init_i_array(n);
    int ip = integer_partition_cumulative_spectrum(n, &v);
    (*spectrum_dest).length = 0;
    append_i_array(spectrum_dest, 1);
    for (i = 1; i < n; i++) {
        append_i_array(spectrum_dest, (v.a[i] - v.a[i-1]));
    }
    free_i_array(&v);
    return ip;
}
    
double integer_partition_cumulative_probs(int n, d_array * probs) {
    int i;
    i_array v;
    v = init_i_array(n);
    int ip = integer_partition_cumulative_spectrum(n, &v);
    (*probs).length = 0;
    for (i = 0; i < n; i++) {
        append_d_array(probs, (v.a[i] / ((double) ip)));
    }
    free_i_array(&v);
    return (get_d_array(probs, (n-1)));
}

double integer_partition_probs(int n, d_array * probs) {
    int i;
    i_array v;
    v = init_i_array(n);
    int ip = integer_partition_spectrum(n, &v);
    (*probs).length = 0;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        append_d_array(probs, (v.a[i] / ((double) ip)));
        sum += get_d_array(probs, i);
    }
    free_i_array(&v);
    return sum;
}

int integer_partition(int n) {
    i_array v;
    v = init_i_array(n);
    int ip = integer_partition_cumulative_spectrum(n, &v);
    free_i_array(&v);
    return ip;
}

int draw_integer_partition_category(const gsl_rng * rng, int n) {
    d_array cumulative_probs;
    cumulative_probs = init_d_array(n);
    double total_prob = integer_partition_cumulative_probs(n,
            &cumulative_probs);
    double r = gsl_rng_uniform(rng);
    int i;
    for (i = 0; i < n; i++) {
        if (r < get_d_array(&cumulative_probs, i)) {
            break;
        }
    }
    free_d_array(&cumulative_probs);
    return (i+1);
}

/** 
 * A function for generating a random draw from a Dirichlet process.
 */
int dirichlet_process_draw(const gsl_rng * rng, int n, double alpha,
        i_array * elements) {
    int num_subsets;
    double subset_prob, new_subset_prob, u;
    i_array subset_counts;
    (*elements).length = 0;
    append_i_array(elements, 0);
    subset_counts = init_i_array(n);
    append_i_array(&subset_counts, 1);
    num_subsets = 1;
    int i, j;
    for (i = 1; i < n; i++) {
        new_subset_prob = (alpha / (alpha + (double)i));
        u = gsl_rng_uniform(rng);
        u -= new_subset_prob;
        if (u < 0.0) {
            append_i_array(elements, num_subsets);
            append_i_array(&subset_counts, 1);
            num_subsets += 1;
            continue;
        }
        for (j = 0; j < num_subsets; j++) {
            subset_prob = ((double)subset_counts.a[j] / (alpha + (double)i));
            u -= subset_prob;
            if (u < 0.0) {
                append_i_array(elements, j);
                subset_counts.a[j] += 1;
                break;
            }
        }
        if (u > 0.0) {
            append_i_array(elements, (num_subsets-1));
            subset_counts.a[num_subsets-1] += 1;
        }
    }
    free_i_array(&subset_counts);
    return num_subsets;
}

/** 
 * A function for generating all partitions of an integer.
 */ 
void generate_integer_partitions(int n, int ip, int ** partitions) {
    /* memset(partitions, 0, sizeof(partitions[0][0]) * ip * n); */
    int i, j;
    for (i = 0; i < ip; i++) {
        for (j = 0; j < n; j++) {
            partitions[i][j] = 0;
        }
    }
    int x[n];
    for (i = 0; i < n; i++) {
        x[i] = 1;
    }
    x[0] = n;
    int m = 0;
    int h = 0;
    int r, t, sum;
    partitions[0][0] = x[0];
    int part_index = 0;
    while (x[0] != 1) {
        part_index += 1;
        if (x[h] == 2) {
            m += 1;
            x[h] = 1;
            h -= 1;
        }
        else {
            r = x[h] - 1;
            t = m - h + 1;
            x[h] = r;
            while (t >= r) {
                h += 1;
                x[h] = r;
                t -= r;
            }
            if (t == 0) {
                m = h;
            }
            else {
                m = h + 1;
                if (t > 1) {
                    h += 1;
                    x[h] = t;
                }
            }
        }
        sum = 0;
        for (i = 0; i < n; i++) {
            partitions[part_index][i] = x[i];
            sum += x[i];
            if (sum >= n) {
                break;
            }
        }
    }
}

