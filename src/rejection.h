/**
 * \file rejection.h
 * \brief A collection of functions for Euclidean distance rejection.
 *
 *
 */

#ifndef REJECTION_H
#define REJECTION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h> // for getopt
#include <string.h>

#define VERSION "0.1"

typedef struct d_array_ {
    double * a;
    int length;
    int capacity;
} d_array;

typedef struct i_array_ {
    int * a;
    int length;
    int capacity;
} i_array;

typedef struct c_array_ {
    char * a;
    int capacity;
} c_array;

typedef struct s_array_ {
    char ** a;
    int length;
    int capacity;
} s_array;

typedef struct sample_sum_ {
    int n;
    double sum;
    double sum_of_squares;
} sample_sum;

typedef struct sample_sum_array_ {
    sample_sum * a;
    int length;
} sample_sum_array;

typedef struct config_ {
    char * observed_path;
    s_array sim_paths;
    int num_retain;
    int num_subsample;
    int means_provided;
    int std_devs_provided;
    d_array means;
    d_array std_devs;
    int include_distance;
} config;

typedef struct sample_ {
    char * file_path;
    int line_num;
    double distance;
    s_array line_array;
} sample;

typedef struct sample_array_ {
    sample * a;
    int length;
    int capacity;
    s_array header;
} sample_array;

d_array init_d_array(int length);
/* void init_d_array(d_array * v, int length); */
void expand_d_array(d_array * v);
void append_d_array(d_array * v, double x);
void extend_d_array(d_array * dest, const d_array * to_add);
double get_d(const d_array * v, int index);
void write_d_array(const d_array * v);
void free_d_array(d_array * v);

c_array init_c_array(int length);
/* void init_c_array(c_array * v, int length); */
void expand_c_array(c_array * v);
void free_c_array(c_array * v);

i_array init_i_array(int length);
/* void init_i_array(i_array * v, int length); */
void expand_i_array(i_array * v);
void append_i_array(i_array * v, int x);
void extend_i_array(i_array * dest, const i_array * to_add);
int get_i(const i_array * v, int index);
void write_i_array(const i_array * v);
void free_i_array(i_array * v);

s_array init_s_array(int length);
/* void init_s_array(s_array * v, int length); */
void expand_s_array(s_array * v);
void append_s_array(s_array * v, const char * x);
void extend_s_array(s_array * dest, const s_array * to_add);
char * get_s(const s_array * v, int index);
void write_s_array(const s_array * v);
void free_s_array(s_array * v);

config init_config();
/* void init_config(config * c); */
void free_config(config * c);
sample init_sample(
        char * file_path,
        const int line_num,
        const s_array * line_array,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs);
/* void init_sample(sample * s, */
/*         char * file_path, */
/*         const int line_num, */
/*         const s_array * line_array, */
/*         const i_array * stat_indices, */
/*         const d_array * std_observed_stats, */
/*         const d_array * means, */
/*         const d_array * std_devs); */
void free_sample(sample * s);
void write_sample(const sample * s, const int include_distance);
sample_array init_sample_array(int length);
/* void init_sample_array(sample_array * v, int length); */
void free_sample_array(sample_array * v);
int process_sample(sample_array * samples, const sample * s);
void rshift_samples(sample_array * s, int index);
void write_sample_array(const sample_array * s, const int include_distance);
sample_sum init_sample_sum();
/* void init_sample_sum(sample_sum * ss); */
sample_sum_array init_sample_sum_array(int length);
/* void init_sample_sum_array(sample_sum_array * s, int length); */
void free_sample_sum_array(sample_sum_array * v);
void update_sample_sum(sample_sum * s, double x);
double get_mean(const sample_sum * s);
double get_sample_variance(const sample_sum * s);
double get_std_dev(const sample_sum * s);
void update_sample_sum_array(sample_sum_array * s, const d_array * x);
void get_mean_array(const sample_sum_array * s, d_array * means);
void get_sample_variance_array(const sample_sum_array * s,
        d_array * v);
void get_std_dev_array(const sample_sum_array * s, d_array * std_devs);
double get_euclidean_distance(const d_array * v1, const d_array * v2);
void standardize_vector(d_array * v, const d_array * means,
        const d_array * std_devs);
void help();
void print_config(const config * c);
void parse_args(config * conf, int argc, char **argv);
int split_str(const c_array * string, s_array * words, int expected_num);
void parse_header(const char * path, c_array * line_buffer, s_array * header);
int headers_match(const s_array * h1, const s_array * h2);
void parse_observed_stats_file(const char * path, c_array * line_buffer,
        s_array * header, d_array * stats);
void get_matching_indices(const s_array * search_strings,
        const s_array * target_strings,
        i_array * indices);
sample_array reject(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        d_array * std_observed_stats,
        d_array * means,
        d_array * std_devs,
        int num_retain,
        s_array * header);
int get_stats(const s_array * line_array, const i_array * stat_indices,
        d_array * stats);
void summarize_stat_samples(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        sample_sum_array * ss_array,
        d_array * means,
        d_array * std_devs,
        int num_to_sample,
        int expected_num_columns);
int main(int argc, char **argv);

#endif

