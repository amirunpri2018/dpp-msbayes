/**
 * @file        eureject.h
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of functions for Euclidean-distance-based
 *              rejection.
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

#ifndef EUREJECT_H
#define EUREJECT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h> // for getopt
#include <string.h>
#include <ctype.h>

#include "math_utils.h"
#include "stats_utils.h"
#include "array_utils.h"
#include "parsing.h"

#define VERSION "0.1"

typedef struct config_ {
    c_array * observed_path;
    s_array * sim_paths;
    int num_retain;
    int num_subsample;
    int means_provided;
    int std_devs_provided;
    d_array * means;
    d_array * std_devs;
    int include_distance;
} config;

typedef struct sample_ {
    c_array * file_path;
    int line_num;
    double distance;
    s_array * line_array;
} sample;

typedef struct sample_array_ {
    sample ** a;
    int length;
    int capacity;
    s_array * header;
} sample_array;

config * init_config();
void free_config(config * c);
sample * init_sample(
        char * file_path,
        const int line_num,
        const s_array * line_array,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs);
void free_sample(sample * s);
void write_sample(const sample * s, const int include_distance);
sample_array * init_sample_array(int length);
void free_sample_array(sample_array * v);
int process_sample(sample_array * samples, sample * s);
void rshift_samples(sample_array * s, int index);
void write_sample_array(const sample_array * s, const int include_distance);
void help();
void print_config(const config * c);
void parse_args(config * conf, int argc, char ** argv);
sample_array * reject(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        d_array * std_observed_stats,
        d_array * means,
        d_array * std_devs,
        int num_retain,
        s_array * header);
void summarize_stat_samples(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        sample_sum_array * ss_array,
        d_array * means,
        d_array * std_devs,
        int num_to_sample,
        int expected_num_columns);
int eureject_main(int argc, char ** argv);

#endif /* EUREJECT_H */

