/**
 * @file        eureject.c
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

#include "eureject.h"

config * init_config() {
    config * c;
    c = (typeof(*c) *) malloc(sizeof(*c));
    c->num_retain = 1000;
    c->num_subsample = 10000;
    c->means_provided = 0;
    c->std_devs_provided = 0;
    c->means = init_d_array(1);
    c->std_devs = init_d_array(1);
    c->sim_paths = init_s_array(1);
    c->observed_path = init_c_array(63);
    c->include_distance = 0;
    return c;
}

void free_config(config * c) {
    free_d_array(c->means);
    free_d_array(c->std_devs);
    free_s_array(c->sim_paths);
    free_c_array(c->observed_path);
    free(c);
    c = NULL;
}

sample * init_sample(
        char * file_path,
        const int line_num,
        const s_array * line_array,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs) {
    int get_stats_return;
    d_array * stats;
    sample * s;
    s = (typeof(*s) *) malloc(sizeof(*s));
    s->file_path = init_c_array(63);
    assign_c_array(s->file_path, file_path);
    s->line_num = line_num;
    s->line_array = init_s_array(line_array->length);
    extend_s_array(s->line_array, line_array);
    stats = init_d_array(stat_indices->length);
    get_stats_return = get_doubles(line_array, stat_indices, stats);
    if (get_stats_return != 0) {
        fprintf(stderr, "ERROR: file %s line %d contains %d invalid stats "
                "columns\n",
                file_path, line_num, get_stats_return);
    }
    standardize_vector(stats, means, std_devs);
    s->distance = get_euclidean_distance(std_observed_stats, stats);
    free_d_array(stats);
    return s;
}

void write_sample(const sample * s, const int include_distance) {
    if (include_distance != 0) {
        fprintf(stdout, "%lf\t", s->distance);
    }
    write_s_array(s->line_array);
}

void free_sample(sample * s) {
    free_s_array(s->line_array);
    free_c_array(s->file_path);
    free(s);
    s = NULL;
}
    
sample_array * init_sample_array(int capacity) {
    if (capacity < 1) {
        fprintf(stderr, "ERROR: init_d_array: capacity must be positive int "
                "greater than 0\n");
        exit(1);
    }
    sample_array * v;
    v = (typeof(*v) *) malloc(sizeof(*v));
    v->capacity = capacity;
    if ((v->a = (typeof(*v->a) *) calloc(v->capacity,
            sizeof(*v->a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v->length = 0;
    return v;
}

int process_sample(sample_array * samples, sample * s) {
    if (samples->length == 0) {
        samples->a[0] = s;
        samples->length++;
        return 0;
    }
    if (samples->a[(samples->length - 1)]->distance <= s->distance) {
        if (samples->length >= samples->capacity) {
            return -1;
        }
        else {
            samples->a[samples->length] = s;
            samples->length++;
            return (samples->length - 1);
        }
    }
    int i;
    for (i = 0; i < samples->length; i++) {
        if (samples->a[i]->distance > s->distance) {
            rshift_samples(samples, i);
            samples->a[i] = s;
            return i;
        }
    }
    return -1;
}

void rshift_samples(sample_array * s, int index) {
    int i, inc;
    inc = 0;
    if (s->length < s->capacity) {
        s->a[s->length] = s->a[(s->length - 1)];
        inc = 1;
    }
    for (i = (s->length - 1); i > index; i--) {
        s->a[i] = s->a[i - 1];
    }
    s->length += inc;
}

void write_sample_array(const sample_array * s, const int include_distance) {
    int i;
    if (include_distance != 0) {
        fprintf(stdout, "distance\t");
    }
    write_s_array(s->header);
    for (i = 0; i < s->length; i++) {
        write_sample(s->a[i], include_distance);
    }
}

void free_sample_array(sample_array * v) {
    int i;
    for (i = 0; i < v->length; i++) {
        free_sample(v->a[i]);
    }
    free(v->a);
    free_s_array(v->header);
    free(v);
    v = NULL;
}

void help() {
    char *version = VERSION;
    printf("EuReject Version %s\n\n", version);
    printf("Usage:\n");
    printf("  eureject -o OBSERVED-FILE [ -k INT -n INT -d ] \\\n");
    printf("      [ -m MEAN1,MEAN2,... -s STDEV1,STDEV2,... ] \\\n");
    printf("      SIMS-FILE1 [ SIMS-FILE2 [...] ]\n\n");
    printf("Options:\n");
    printf(" -f  Path to file containing observed summary statistics\n");
    printf(" -k  Number of samples to keep. Default: 1000.\n");
    printf("     If set to 0, only the stat means and standard deviations\n");
    printf("     are calculated and reported (i.e., no rejection is\n");
    printf("     performed).\n");
    printf(" -n  Number of samples to use for calculating stat means and\n");
    printf("     standard deviations for standardizing the statistics.\n");
    printf("     This option is ignored if `-m` and `-s` are provided.\n");
    printf("     Default: 10000\n");
    printf(" -m  Comma-separated list of means to use for standardizing\n");
    printf("     statistics. The number of means must match the number of\n");
    printf("     columns in the observed stats file, and must be in the\n");
    printf("     same order. Must be used with `-s`.\n");
    printf(" -s  Comma-separated list of standard deviations to use for\n");
    printf("     standardizing statistics. The number of means must match\n");
    printf("     the number of columns in the observed stats file, and must\n");
    printf("     be in the same order. Must be used with `-m`.\n");
    printf(" -e  Report Euclidean distances of retained samples in the\n");
    printf("     first column of the output. Default is not to report\n");
    printf("     the distance column.\n");
    printf(" -h  Display this help message and exit\n");
}

void print_config(const config * c) {
    int i;
    fprintf(stderr, "Number of samples to retain: %d\n", c->num_retain);
    fprintf(stderr, "Number of samples to use for standardization: %d\n",
            c->num_subsample);
    fprintf(stderr, "Observed stats path: %s\n", c->observed_path->a);
    fprintf(stderr, "Path(s) to file(s) with simulated draws: ");
    for (i = 0; i < (c->sim_paths->length - 1); i++) {
        fprintf(stderr, "%s, ", get_s_array(c->sim_paths, i));
    }
    fprintf(stderr, "%s\n", get_s_array(c->sim_paths,
                (c->sim_paths->length-1)));
    fprintf(stderr, "Means for standardization: ");
    if (c->means_provided == 0) {
        fprintf(stderr, "None\n");
    }
    else {
        for (i = 0; i < (c->means->length - 1); i++) {
            fprintf(stderr, "%f, ", c->means->a[i]);
        }
        fprintf(stderr, "%f\n", c->means->a[(c->means->length-1)]);
    }
    fprintf(stderr, "Standard deviations for standardization: ");
    if (c->std_devs_provided == 0) {
        fprintf(stderr, "None\n");
    }
    else {
        for (i = 0; i < (c->std_devs->length - 1); i++) {
            fprintf(stderr, "%f, ", c->std_devs->a[i]);
        }
        fprintf(stderr, "%f\n", c->std_devs->a[(c->std_devs->length-1)]);
    }
}

void parse_args(config * conf, int argc, char ** argv) {
    int i, j;
    char * p;
    char * end_ptr;
    char * end_ptr_orig;
    /* opterr = 0; */
    end_ptr = (typeof(*end_ptr) *) malloc(sizeof(end_ptr) * 64);
    end_ptr_orig = end_ptr;
    conf->means->length = 0;
    conf->std_devs->length = 0;
    conf->sim_paths->length = 0;
    while((i = getopt(argc, argv, "f:k:n:m:s:eh")) != -1) {
        switch(i) {
            case 'f':
                assign_c_array(conf->observed_path, optarg);
                break;
            case 'k':
                conf->num_retain = atoi(optarg);
                break;
            case 'n':
                if (atoi(optarg) > 0) {
                    conf->num_subsample = atoi(optarg);
                    break;
                }
                else {
                    fprintf(stderr, "ERROR: `-n' must be positive integer\n");
                    help();
                    exit(1);
                }
                break;
            case 'm':
                j = 0;
                conf->means_provided = 1;
                p = strtok(optarg, ",");
                while (p != NULL) {
                    append_d_array(conf->means, strtod(p, &end_ptr));
                    if (end_ptr == p) {
                        fprintf(stderr, "ERROR: mean %d is not a valid "
                                "number\n", (j + 1));
                    }
                    p = strtok(NULL, ",");
                    j += 1;
                }
                break;
            case 's':
                j = 0;
                conf->std_devs_provided = 1;
                p = strtok(optarg, ",");
                while (p != NULL) {
                    append_d_array(conf->std_devs, strtod(p, &end_ptr));
                    if (end_ptr == p) {
                        fprintf(stderr, "ERROR: standard deviation %d is not "
                                "a valid number\n", (j + 1));
                    }
                    p = strtok(NULL, ",");
                    j += 1;
                }
                break;
            case 'e':
                conf->include_distance = 1;
                break;
            case 'h':
                help();
                exit(0);
                break;
            case '?':
                if (optopt == 'o') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'k') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'n') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 'm') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (optopt == 's') {
                    fprintf(stderr, "ERROR: option `-%c' requires an "
                            "argument\n", optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "ERROR: unknown option `-%c'\n", optopt);
                }
                else {
                    fprintf(stderr, "ERROR: unknown option character `\\x%x'\n",
                            optopt);
                }
                break;
            default:
                help();
                exit(1);
        }
    }
    for (i = optind, j = 0; i < argc; i++, j++) {
        append_s_array(conf->sim_paths, argv[i]);
    }
    if (conf->means_provided == 1) {
        conf->num_subsample = 0;
    }
    // vetting
    if (conf->observed_path->a == NULL) {
        fprintf(stderr, "ERROR: Please provide path to observed stats\n");
        help();
        exit(1);
    }
    if ((conf->sim_paths->length < 1) ||
            (get_s_array(conf->sim_paths, 0) == NULL)) {
        fprintf(stderr, "ERROR: Please provide at least one simulation file\n");
        help();
        exit(1);
    }
    if (((conf->means_provided == 0) && (conf->std_devs_provided == 1)) ||
            ((conf->means_provided == 1) && (conf->std_devs_provided == 0))) {
        fprintf(stderr, "ERROR: Please specify both means and standard "
                "deviations for standardization, or neither.\n");
        help();
        exit(1);
    }
    if (conf->means_provided && (conf->means->length != conf->std_devs->length)) {
        fprintf(stderr, "ERROR: Please provide equal numbers of means and "
                "std deviations\n");
        help();
        exit(1);
    }
    free(end_ptr_orig);
}

sample_array * reject(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        d_array * std_observed_stats,
        d_array * means,
        d_array * std_devs,
        int num_retain,
        s_array * header) {
    FILE * f;
    int i, line_num, ncols, sample_idx;
    s_array * line_array;
    sample_array * retained_samples;
    line_array = init_s_array((*header).length);
    retained_samples = init_sample_array(num_retain);
    retained_samples->header = init_s_array((*header).length);
    extend_s_array(retained_samples->header, header);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen(get_s_array(paths, i), "r")) == NULL) {
            perror(get_s_array(paths, i));
            exit(1);
        }
        while (fgets((*line_buffer).a, (((*line_buffer).capacity) - 1),
                    f) != NULL) {
            line_num++;
            ncols = split_str((*line_buffer).a, line_array, (*header).length);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", get_s_array(paths, i), line_num,
                        ncols, (*header).length);
                exit(1);
            }
            if (line_num == 1) continue;
            sample * s;
            s = init_sample(get_s_array(paths, i), line_num, line_array,
                    stat_indices, std_observed_stats, means, std_devs);
            sample_idx = process_sample(retained_samples, s);
            if (sample_idx < 0) {
                free_sample(s);
            }
        }
        fclose(f);
    }
    free_s_array(line_array);
    return retained_samples;
}

void summarize_stat_samples(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        sample_sum_array * ss_array,
        d_array * means,
        d_array * std_devs,
        int num_to_sample,
        int expected_num_columns) {
    FILE * f;
    int i, line_num, ncols, get_stats_return;
    s_array * line_array;
    d_array * stats;
    line_array = init_s_array(expected_num_columns);
    stats = init_d_array((*stat_indices).length);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen(get_s_array(paths, i), "r")) == NULL) {
            perror(get_s_array(paths, i));
            exit(1);
        }
        while (fgets((*line_buffer).a, (((*line_buffer).capacity) - 1),
                    f) != NULL) {
            line_num++;
            ncols = split_str((*line_buffer).a, line_array,
                    expected_num_columns);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", get_s_array(paths, i), line_num,
                        ncols, expected_num_columns);
                exit(1);
            }
            if (line_num == 1) continue;
            get_stats_return = get_doubles(line_array, stat_indices, stats);
            if (get_stats_return != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d invalid stats "
                        "columns\n", get_s_array(paths, i), line_num,
                        get_stats_return);
                exit(1);
            }
            update_sample_sum_array(ss_array, stats);
            if (ss_array->a[0]->n >= num_to_sample) break;
        }
        fclose(f);
        if (ss_array->a[0]->n >= num_to_sample) break;
    }
    get_mean_array(ss_array, means);
    get_std_dev_array(ss_array, std_devs);
    free_d_array(stats);
    free_s_array(line_array);
}

int eureject_main(int argc, char ** argv) {
    c_array * line_buffer;
    s_array * obs_header;
    s_array * sim_header;
    s_array * sim_header_comp;
    d_array * obs_stats;
    int i, heads_match;
    i_array * indices;
    sample_sum_array * sample_sums;
    sample_array * retained_samples;
    config * conf;
    line_buffer = init_c_array(pow(2, 10));
    obs_header = init_s_array(1);
    obs_stats = init_d_array(1);
    if (argc < 2) {
        help();
        exit(1);
    }
    conf = init_config();
    parse_args(conf, argc, argv);
    print_config(conf);

    parse_observed_stats_file(conf->observed_path->a, line_buffer, obs_header,
            obs_stats);
    sim_header = init_s_array(obs_header->length);
    parse_header(get_s_array(conf->sim_paths, 0), line_buffer, sim_header);
    if (conf->sim_paths->length > 1) {
        sim_header_comp = init_s_array(sim_header->length);
        for (i = 1; i < conf->sim_paths->length; i++) {
            parse_header(get_s_array(conf->sim_paths, i), line_buffer,
                    sim_header_comp);
            heads_match = s_arrays_equal(sim_header, sim_header_comp);
            if (heads_match != 1) {
                fprintf(stderr, "ERROR: Files %s and %s have different "
                        "headers\n", get_s_array(conf->sim_paths, 0),
                        get_s_array(conf->sim_paths, i));
                exit(1);
            }
        }
        free_s_array(sim_header_comp);
    }
    indices = init_i_array(obs_header->length);
    get_matching_indices(obs_header, sim_header, indices);
    if (conf->means_provided == 0) {
        sample_sums = init_sample_sum_array(obs_header->length);
        summarize_stat_samples(conf->sim_paths, line_buffer, indices,
                sample_sums, conf->means, conf->std_devs, conf->num_subsample,
                sim_header->length);
        free_sample_sum_array(sample_sums);
    }
    if (conf->num_retain < 1) {
        write_s_array(obs_header);
        write_d_array(conf->means);
        write_d_array(conf->std_devs);
        free_i_array(indices);
        free_c_array(line_buffer);
        free_s_array(obs_header);
        free_s_array(sim_header);
        free_d_array(obs_stats);
        free_config(conf);
        return 0;
    }
    standardize_vector(obs_stats, conf->means, conf->std_devs);
    retained_samples = reject(conf->sim_paths, line_buffer, indices,
            obs_stats, conf->means, conf->std_devs, conf->num_retain,
            sim_header);
    write_sample_array(retained_samples, conf->include_distance);
    free_sample_array(retained_samples);
    free_i_array(indices);
    free_c_array(line_buffer);
    free_s_array(obs_header);
    free_s_array(sim_header);
    free_d_array(obs_stats);
    free_config(conf);
    return 0;
}
