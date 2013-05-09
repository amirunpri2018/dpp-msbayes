/**
 * \file rejection.c
 * \brief A collection of functions for Euclidean distance rejection.
 *
 *
 */

#include "rejection.h"

sample_sum init_sample_sum() {
    sample_sum ss;
    ss.n = 0;
    ss.sum = 0.0;
    ss.sum_of_squares = 0.0;
    return ss;
}

d_array init_d_array(int length) {
    d_array v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc(v.length, sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    return v;
}

void free_d_array(d_array * v) {
    free((*v).a);
}

void append_d_array(d_array * v, double x) {
    if (((*v).a = (typeof(*(*v).a) *) realloc((*v).a ,
            (((*v).length + 1) * sizeof(*(*v).a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    (*v).a[(*v).length] = x;
    (*v).length += 1;
}

i_array init_i_array(int length) {
    i_array v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc(v.length, sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    return v;
}

void free_i_array(i_array * v) {
    free((*v).a);
}

void append_i_array(i_array * v, int x) {
    if (((*v).a = (typeof(*(*v).a) *) realloc((*v).a ,
            (((*v).length + 1) * sizeof(*(*v).a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    (*v).a[(*v).length] = x;
    (*v).length += 1;
}

c_array init_c_array(int length) {
    c_array v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc((v.length + 1), sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    v.a[v.length] = '\0';
    return v;
}

void free_c_array(c_array * v) {
    free((*v).a);
}

void append_c_array(c_array * v, char x) {
    if (((*v).a = (typeof(*(*v).a) *) realloc((*v).a ,
            (((*v).length + 2) * sizeof(*(*v).a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    (*v).a[(*v).length] = x;
    (*v).a[((*v).length + 1)] = '\0';
    (*v).length += 1;
}

c_array_2d init_c_array_2d(int width, int length) {
    c_array_2d v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc((v.length), sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    int i;
    for (i = 0; i < v.length; i++) {
        v.a[i] = init_c_array(width);
    }
    return v;
}

void free_c_array_2d(c_array_2d * v) {
    int i;
    for (i = 0; i < (*v).length; i++) {
        free(&(*v).a[i]);
    }
    free((*v).a);
}

void append_c_array_2d(c_array_2d * v, c_array * x) {
    if (((*v).a = (typeof(*(*v).a) *) realloc((*v).a ,
            (((*v).length + 1) * sizeof(*(*v).a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    (*v).a[(*v).length] = *x;
    (*v).length += 1;
}

s_array init_s_array(int length) {
    s_array v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc(v.length, sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    return v;
}

void free_s_array(s_array * v) {
    int i;
    for (i = 0; i < (*v).length; i++) {
        free((*v).a[i]);
    }
    free((*v).a);
}

void append_s_array(s_array * v, char * x) {
    if (((*v).a = (typeof(*(*v).a) *) realloc((*v).a ,
            (((*v).length + 1) * sizeof(*(*v).a)))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    (*v).a[(*v).length] = x;
    (*v).length += 1;
}

config init_config() {
    config c;
    c.num_retain = 1000;
    c.num_subsample = 10000;
    c.means_provided = 0;
    c.std_devs_provided = 0;
    c.means = init_d_array(1);
    c.std_devs = init_d_array(1);
    c.sim_paths = init_s_array(1);
    return c;
}

void free_config(config * c) {
    free_d_array(&(*c).means);
    free_d_array(&(*c).std_devs);
    free_s_array(&(*c).sim_paths);
}

sample init_sample(
        char * file_path,
        const int line_num,
        const s_array * line_array,
        const i_array * stat_indices,
        const d_array * std_observed_stats,
        const d_array * means,
        const d_array * std_devs) {
    int i, get_stats_return;
    d_array stats;
    sample s;
    s.file_path = file_path;
    s.line_num = line_num;
    s.line_array = *line_array;
    stats = init_d_array((*stat_indices).length);
    get_stats_return = get_stats(line_array, stat_indices, &stats);
    if (get_stats_return != 0) {
        fprintf(stderr, "ERROR: file %s line %d contains %d invalid stats "
                "columns\n",
                file_path, line_num, get_stats_return);
    }
    standardize_vector(&stats, means, std_devs);
    s.distance = get_euclidean_distance(std_observed_stats, &stats);
    return s;
}

void free_sample(sample * s) {
    int i;
    free_s_array(&(*s).line_array);
}
    
sample_array init_sample_array(int length) {
    sample_array v;
    v.length = length;
    v.sample_size = 0;
    if ((v.a = (typeof(*v.a) *) calloc(v.length, sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    return v;
}

void free_sample_array(sample_array * v) {
    int i;
    for (i = 0; i < (*v).length; i++) {
        free_sample(&(*v).a[i]);
    }
    free((*v).a);
}

int process_sample(sample_array * samples, const sample * s) {
    if ((*samples).sample_size == 0) {
        (*samples).a[0] = *s;
        (*samples).sample_size++;
        return 0;
    }
    if ((*samples).a[((*samples).sample_size - 1)].distance <= (*s).distance) {
        if ((*samples).sample_size >= (*samples).length) {
            return -1;
        }
        else {
            (*samples).a[(*samples).sample_size] = *s;
            (*samples).sample_size++;
            return ((*samples).sample_size + 1);
        }
    }
    int i;
    for (i = 0; i < (*samples).sample_size; i++) {
        if ((*samples).a[i].distance > (*s).distance) {
            rshift_samples(samples, i);
            (*samples).a[i] = *s;
            return i;
        }
    }
    return -1;
}

void rshift_samples(sample_array * s, int index) {
    int i, inc;
    inc = 0;
    if ((*s).sample_size < (*s).length) {
        (*s).a[(*s).sample_size] = (*s).a[((*s).sample_size - 1)];
        inc = 1;
    }
    for (i = ((*s).sample_size - 1); i > index; i--) {
        (*s).a[i] = (*s).a[i - 1];
    }
    (*s).sample_size += inc;
}


sample_sum_array init_sample_sum_array(int length) {
    sample_sum_array v;
    v.length = length;
    if ((v.a = (typeof(*v.a) *) calloc(v.length, sizeof(*v.a))) == NULL) {
        perror("out of memory");
        exit(1);
    }
    int i;
    for (i = 0; i < v.length; i++) {
        v.a[i] = init_sample_sum();
    }
    return v;
}

void free_sample_sum_array(sample_sum_array * v) {
    free((*v).a);
}

void update_sample_sum(sample_sum * s, double x) {
    (*s).n += 1;
    (*s).sum += x;
    (*s).sum_of_squares += pow(x, 2);
}

double get_mean(const sample_sum * s) {
    if ((*s).n < 1) {
        return 0.0;
    }
    double mn = ((*s).sum / (*s).n);
    return mn;
}

double get_sample_variance(const sample_sum * s) {
    if ((*s).n < 2) {
        return -1.0;
    }
    double mn = get_mean(s);
    double v = (((*s).sum_of_squares - (mn * (*s).sum)) / ((*s).n - 1));
    return v;
}

double get_std_dev(const sample_sum * s) {
    if ((*s).n < 2) {
        return -1.0;
    }
    double var = get_sample_variance(s);
    double sd = (sqrt(var));
    return sd;
}

void update_sample_sum_array(sample_sum_array * s, const d_array * x) {
    if ((*s).length != (*x).length) {
        fprintf(stderr, "ERROR: update_sample_sum_array: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < (*s).length; i++) {
        update_sample_sum(&(*s).a[i], (*x).a[i]);
    }
}

void get_mean_array(const sample_sum_array * s, d_array * means) {
    if ((*s).length != (*means).length) {
        fprintf(stderr, "ERROR: get_mean_array: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < (*s).length; i++) {
        (*means).a[i] = get_mean(&(*s).a[i]);
    }
}

void get_sample_variance_array(const sample_sum_array * s,
        d_array * v) {
    if ((*s).length != (*v).length) {
        fprintf(stderr, "ERROR: get_sample_variance_array: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < (*s).length; i++) {
        (*v).a[i] = get_sample_variance(&(*s).a[i]);
    }
}
        
void get_std_dev_array(const sample_sum_array * s, d_array * std_devs) {
    if ((*s).length != (*std_devs).length) {
        fprintf(stderr, "ERROR: get_sample_variance_array: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < (*s).length; i++) {
        (*std_devs).a[i] = get_std_dev(&(*s).a[i]);
    }
}

double get_euclidean_distance(const d_array * v1, const d_array * v2) {
    if ((*v1).length != (*v2).length) {
        fprintf(stderr, "ERROR: get_euclidean_distance: arrays must be of "
                "equal length\n");
        exit(1);
    }
    double sum_of_squared_diffs = 0.0;
    int i;
    for (i = 0; i < (*v1).length; i++) {
        sum_of_squared_diffs += pow(((*v1).a[i] - (*v2).a[i]), 2);
    }
    return sqrt(sum_of_squared_diffs);
}

void standardize_vector(d_array * v, const d_array * means,
        const d_array * std_devs) {
    if (((*v).length != (*means).length) ||
            ((*v).length != (*std_devs).length)) {
        fprintf(stderr, "ERROR: standardize_vector: arrays must be of "
                "equal length\n");
        exit(1);
    }
    int i;
    for (i = 0; i < (*v).length; i++) {
        (*v).a[i] = ((*v).a[i] - (*means).a[i]) / (*std_devs).a[i];
    }
}

void help() {
    char *version = VERSION;
    printf("EuReject Version %s\n\n", version);
    printf("Usage:\n");
    printf("  eureject -o OBSERVED-FILE [ -k INT -n INT ] \\\n");
    printf("      [ -m MEAN1,MEAN2,... -s STDEV1,STDEV2,... ] \\\n");
    printf("      SIMS-FILE1 [ SIMS-FILE2 [...] ]\n\n");
    printf("Options:\n");
    printf(" -o  Path to file containing observed summary statistics\n");
    printf(" -k  Number of samples to keep. Default: 1000\n");
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
}

void print_config(const config * c) {
    int i;
    fprintf(stderr, "Number of samples to retain: %d\n", (*c).num_retain);
    fprintf(stderr, "Number of samples to use for standardization: %d\n",
            (*c).num_subsample);
    fprintf(stderr, "Observed stats path: %s\n", (*c).observed_path);
    fprintf(stderr, "Path(s) to files with simulated draws: ");
    for (i = 0; i < ((*c).sim_paths.length - 1); i++) {
        fprintf(stderr, "%s, ", (*c).sim_paths.a[i]);
    }
    fprintf(stderr, "%s\n", (*c).sim_paths.a[((*c).sim_paths.length-1)]);
    fprintf(stderr, "Means for standardization: ");
    if ((*c).means_provided == 0) {
        fprintf(stderr, "None\n");
    }
    else {
        for (i = 0; i < ((*c).means.length - 1); i++) {
            fprintf(stderr, "%f, ", (*c).means.a[i]);
        }
        fprintf(stderr, "%f\n", (*c).means.a[((*c).means.length-1)]);
    }
    fprintf(stderr, "Standard deviations for standardization: ");
    if ((*c).std_devs_provided == 0) {
        fprintf(stderr, "None\n");
    }
    else {
        for (i = 0; i < ((*c).std_devs.length - 1); i++) {
            fprintf(stderr, "%f, ", (*c).std_devs.a[i]);
        }
        fprintf(stderr, "%f\n", (*c).std_devs.a[((*c).std_devs.length-1)]);
    }
}

void parse_args(config * conf, int argc, char **argv) {
    int i, j, k;
    char * p;
    /* opterr = 0; */
    while((i = getopt(argc, argv, "o:k:n:m:s:h")) != -1) {
        switch(i) {
            case 'o':
                (*conf).observed_path = optarg;
                break;
            case 'k':
                if (atoi(optarg) > 0) {
                    (*conf).num_retain = atoi(optarg);
                    break;
                }
                else {
                    error("`-k' must be positive integer");
                    help();
                    exit(1);
                }
            case 'n':
                if (atoi(optarg) > 0) {
                    (*conf).num_subsample = atoi(optarg);
                    break;
                }
                else {
                    error("`-n' must be positive integer");
                    help();
                    exit(1);
                }
            case 'm':
                j = 0;
                (*conf).means_provided = 1;
                p = strtok(optarg, ",");
                while (p != NULL) {
                    if (j < (*conf).means.length) {
                        (*conf).means.a[j] = atof(p);
                    }
                    else {
                        append_d_array(&(*conf).means, atof(p));
                    }
                    p = strtok(NULL, ",");
                    j += 1;
                }
                break;
            case 's':
                j = 0;
                (*conf).std_devs_provided = 1;
                p = strtok(optarg, ",");
                while (p != NULL) {
                    if (j < (*conf).std_devs.length) {
                        (*conf).std_devs.a[j] = atof(p);
                    }
                    else {
                        append_d_array(&(*conf).std_devs, atof(p));
                    }
                    p = strtok(NULL, ",");
                    j += 1;
                }
                break;
            case 'h':
                help();
                exit(0);
            case '?':
                if (optopt == 'o') {
                    fprintf(stderr, "option `-%c' requires an argument\n",
                            optopt);
                }
                else if (optopt == 'k') {
                    fprintf(stderr, "option `-%c' requires an argument\n",
                            optopt);
                }
                else if (optopt == 'n') {
                    fprintf(stderr, "option `-%c' requires an argument\n",
                            optopt);
                }
                else if (optopt == 'm') {
                    fprintf(stderr, "option `-%c' requires an argument\n",
                            optopt);
                }
                else if (optopt == 's') {
                    fprintf(stderr, "option `-%c' requires an argument\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "unknown option `-%c'\n", optopt);
                }
                else {
                    fprintf(stderr, "unknown option character `\\x%x'\n",
                            optopt);
                }
            default:
                help();
                exit(1);
        }
    }
    for (i = optind, j = 0; i < argc; i++, j++) {
        if (j < (*conf).sim_paths.length) {
            (*conf).sim_paths.a[j] = argv[i];
        }
        else {
            append_s_array(&(*conf).sim_paths, argv[i]);
        }
    }
}

int split_str(c_array * string, s_array * words, int expected_num) {
    int column_idx, n;
    char * ptr;
    c_array match;
    column_idx = 0;
    ptr = (*string).a;
    while(*ptr) {
        match = init_c_array(63);
        if ((sscanf(ptr, "%s%n", match.a, &n)) == 1) {
            ptr += n;
            if (column_idx < (*words).length) {
                (*words).a[column_idx] = match.a;
            }
            else {
                append_s_array(words, match.a);
            }
            column_idx++;
        }
        ++ptr;
    }
    if ((expected_num > 0) && (expected_num != column_idx)) {
        if (column_idx == 0) column_idx--;
        return column_idx;
    }
    return 0;
    // free match
}

void parse_header(const char * path, c_array * line_buffer, s_array * header) {
    FILE * f;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    // parse header
    if ((fgets((*line_buffer).a, (((*line_buffer).length) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no lines in %s", path);
        exit(1);
    }
    split_str(line_buffer, header, 0);
    fclose(f);
}

int headers_match(const s_array * h1, const s_array * h2) {
    int i;
    if ((*h1).length != (*h2).length) {
        return 0;
    }
    for (i = 0; i < (*h1).length; i++) {
        if (strcmp((*h1).a[i], (*h2).a[i]) != 0) {
            return 0;
        }
    }
    return 1;
}

void parse_observed_stats_file(const char * path, c_array * line_buffer,
        s_array * header, d_array * stats) {
    FILE * f;
    int column_idx, n;
    char * ptr;
    double fmatch;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    // parse header
    if ((fgets((*line_buffer).a, (((*line_buffer).length) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no header in %s", path);
        exit(1);
    }
    split_str(line_buffer, header, 0);
    // parse stats
    if ((fgets((*line_buffer).a, (((*line_buffer).length) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no stats in %s", path);
        exit(1);
    }
    column_idx = 0;
    ptr = (*line_buffer).a;
    while(*ptr) {
        if ((sscanf(ptr, "%lf%n", &fmatch, &n)) == 1) {
            ptr += n;
            if (column_idx < (*stats).length) {
                (*stats).a[column_idx] = fmatch;
            }
            else {
                append_d_array(stats, fmatch);
            }
            column_idx++;
        }
        ++ptr;
    }
    fclose(f);
}

void get_matching_indices(const s_array * search_strings,
        const s_array * target_strings,
        i_array * indices) {
    int i, j, found;
    for (i = 0; i < (*search_strings).length; i++) {
        found = 0;
        for (j = 0; j < (*target_strings).length; j++) {
            if (strcmp((*search_strings).a[i], (*target_strings).a[j]) == 0) {
                if (found != 0) {
                    fprintf(stderr, "string %s found more than once\n",
                            (*search_strings).a[i]);
                    exit(1);
                }
                found = 1;
                if (i < (*indices).length) {
                    (*indices).a[i] = j;
                }
                else {
                    append_i_array(indices, j);
                }
            }
        }
        if (found == 0) {
            fprintf(stderr, "string %s was not found\n",
                    (*search_strings).a[i]);
            exit(1);
        }
    }
}

sample_array reject(const s_array * paths,
        c_array * line_buffer,
        const i_array * stat_indices,
        d_array * std_observed_stats,
        d_array * means,
        d_array * std_devs,
        int num_retain,
        int expected_num_columns) {
    FILE * f;
    int i, j, line_num, ncols, get_stats_return, sample_idx;
    s_array line_array;
    sample_array retained_samples;
    retained_samples = init_sample_array(num_retain);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen((*paths).a[i], "r")) == NULL) {
            perror((*paths).a[i]);
            exit(1);
        }
        while (fgets((*line_buffer).a, (((*line_buffer).length) - 1),
                    f) != NULL) {
            line_num++;
            line_array = init_s_array(expected_num_columns);
            ncols = split_str(line_buffer, &line_array, expected_num_columns);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", (*paths).a[i], line_num, ncols,
                        expected_num_columns);
                exit(1);
            }
            if (line_num == 1) continue;
            sample s;
            s = init_sample((*paths).a[i], line_num, &line_array, stat_indices,
                    std_observed_stats, means, std_devs);
            sample_idx = process_sample(&retained_samples, &s);
            fprintf(stderr, "file %s line %d: %lf\n", s.file_path, s.line_num,
                    s.distance);
            if (sample_idx < 0) {
                free_sample(&s);
            }
        }
        fclose(f);
    }
    return retained_samples;
}

int get_stats(const s_array * line_array, const i_array * stat_indices,
        d_array * stats) {
    int i, ret;
    char * end_ptr;
    end_ptr = (typeof(*end_ptr) *) malloc(sizeof(end_ptr) * 64);
    ret = 0;
    for (i = 0; i < (*stat_indices).length; i++) {
        (*stats).a[i] = strtod((*line_array).a[(*stat_indices).a[i]], &end_ptr);
        if (end_ptr == (*line_array).a[(*stat_indices).a[i]]) {
            fprintf(stderr, "ERROR: column %d is not a valid "
                    "number\n", ((*stat_indices).a[i] + 1));
            ret++;
        }
    }
    return ret;
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
    int i, j, line_num, ncols, get_stats_return;
    s_array line_array;
    d_array stats;
    line_array = init_s_array(expected_num_columns);
    stats = init_d_array((*stat_indices).length);
    for (i = 0; i < (*paths).length; i++) {
        line_num = 0;
        if ((f = fopen((*paths).a[i], "r")) == NULL) {
            perror((*paths).a[i]);
            exit(1);
        }
        while (fgets((*line_buffer).a, (((*line_buffer).length) - 1),
                    f) != NULL) {
            line_num++;
            ncols = split_str(line_buffer, &line_array, expected_num_columns);
            if (ncols == -1) continue; //empty line
            if (ncols != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d columns "
                        "(expected %d)\n", (*paths).a[i], line_num, ncols,
                        expected_num_columns);
                exit(1);
            }
            if (line_num == 1) continue;
            get_stats_return = get_stats(&line_array, stat_indices, &stats);
            if (get_stats_return != 0) {
                fprintf(stderr, "ERROR: file %s line %d has %d invalid stats "
                        "columns\n", (*paths).a[i], line_num, get_stats_return);
                exit(1);
            }
            update_sample_sum_array(ss_array, &stats);
            if ((*ss_array).a[0].n >= num_to_sample) break;
        }
        fclose(f);
        if ((*ss_array).a[0].n >= num_to_sample) break;
    }
    get_mean_array(ss_array, means);
    get_std_dev_array(ss_array, std_devs);
    free_d_array(&stats);
    free_s_array(&line_array);
}

int main(int argc, char **argv) {
    c_array line_buffer;
    s_array obs_header, sim_header, sim_header_comp;
    d_array obs_stats;
    int i, heads_match;
    i_array indices;
    sample_sum_array sample_sums;
    sample_array retained_samples;
    line_buffer = init_c_array(pow(2, 10));
    obs_header = init_s_array(1);
    obs_stats = init_d_array(1);
    if (argc < 2) {
        help();
        exit(1);
    }
    for (i = 0; i < argc; i++) {
        printf("%s\n", argv[i]);
    }
    config conf;
    conf = init_config();
    parse_args(&conf, argc, argv);
    print_config(&conf);
    if (conf.observed_path == NULL) {
        fprintf(stderr, "ERROR: Please provide path to observed stats\n");
        help();
        exit(1);
    }
    if ((conf.sim_paths.length < 1) || (conf.sim_paths.a[0] == NULL)) {
        fprintf(stderr, "ERROR: Please provide at least one simulation file\n");
        help();
        exit(1);
    }
    if (((conf.means_provided == 0) && (conf.std_devs_provided == 1)) ||
            ((conf.means_provided == 1) && (conf.std_devs_provided == 0))) {
        fprintf(stderr, "ERROR: Please specify both means and standard "
                "deviations for standardization, or neither.\n");
        help();
        exit(1);
    }
    if (conf.means_provided && (conf.means.length != conf.std_devs.length)) {
        fprintf(stderr, "ERROR: Please provide equal numbers of means and "
                "std deviations\n");
        help();
        exit(1);
    }
    if (conf.means_provided == 1) {
        conf.num_subsample = 0;
    }
    print_config(&conf);

    parse_observed_stats_file(conf.observed_path, &line_buffer, &obs_header,
            &obs_stats);
    for (i = 0; i < obs_header.length; i++) {
        printf("%s\n", obs_header.a[i]);
    }
    for (i = 0; i < obs_stats.length; i++) {
        printf("%lf\n", obs_stats.a[i]);
    }
    sim_header = init_s_array(obs_header.length);
    parse_header(conf.sim_paths.a[0], &line_buffer, &sim_header);
    for (i = 0; i < sim_header.length; i++) {
        printf("%s\n", sim_header.a[i]);
    }
    if (conf.sim_paths.length > 1) {
        sim_header_comp = init_s_array(sim_header.length);
        for (i = 1; i < conf.sim_paths.length; i++) {
            parse_header(conf.sim_paths.a[i], &line_buffer, &sim_header_comp);
            heads_match = headers_match(&sim_header, &sim_header_comp);
            if (heads_match != 1) {
                fprintf(stderr, "ERROR: Files %s and %s have different "
                        "headers\n", conf.sim_paths.a[0],
                        conf.sim_paths.a[i]);
                exit(1);
            }
        }
    }
    indices = init_i_array(obs_header.length);
    get_matching_indices(&obs_header, &sim_header, &indices);
    for (i = 0; i < indices.length; i++) {
        printf("%d\n", indices.a[i]);
    }
    if (conf.means_provided == 0) {
        sample_sums = init_sample_sum_array(obs_header.length);
        conf.means = init_d_array(obs_header.length);
        conf.std_devs = init_d_array(obs_header.length);
        summarize_stat_samples(&conf.sim_paths, &line_buffer, &indices,
                &sample_sums, &conf.means, &conf.std_devs, conf.num_subsample,
                sim_header.length);
    }
    for (i = 0; i < conf.means.length; i++) {
        printf("%lf %lf\n", conf.means.a[i], conf.std_devs.a[i]);
    }
    standardize_vector(&obs_stats, &conf.means, &conf.std_devs);
    retained_samples = reject(&conf.sim_paths, &line_buffer, &indices,
            &obs_stats, &conf.means, &conf.std_devs, conf.num_retain,
            sim_header.length);
    printf("%d %d\n", retained_samples.length, retained_samples.sample_size);
    for (i = 0; i < retained_samples.sample_size; i++) {
        printf("%lf\n", retained_samples.a[i].distance);
    }
    free_sample_array(&retained_samples);
    return 0;
}


    /* int i, l = 4; */
    /* d_array v1, v2, means, std_devs; */
    /* v1 = init_d_array(l); */
    /* v2 = init_d_array(l); */
    /* means = init_d_array(l); */
    /* std_devs = init_d_array(l); */
    /* for (i = 0; i < l; i++) { */
    /*     v1.a[i] = 2; */
    /*     v2.a[i] = 4; */
    /*     means.a[i] = 3; */
    /*     std_devs.a[i] = 2; */
    /* } */
    /* double d1 = get_euclidean_distance(&v1, &v2); */
    /* standardize_vector(&v1, &means, &std_devs); */
    /* standardize_vector(&v2, &means, &std_devs); */
    /* double d2 = get_euclidean_distance(&v1, &v2); */
    /* fprintf(stdout, "%f\n%f\n", d1, d2); */
    /* /1* d_array v3 = init_d_array(5); *1/ */
    /* /1* double d3 = get_euclidean_distance(v1, v3); *1/ */

    /* sample_sum s; */
    /* s = init_sample_sum(); */
    /* update_sample_sum(&s, 2.0); */
    /* update_sample_sum(&s, 4.0); */
    /* update_sample_sum(&s, 6.0); */
    /* double m, v, sd; */
    /* m = get_mean(&s); */
    /* v = get_sample_variance(&s); */
    /* sd = get_std_dev(&s); */
    /* fprintf(stdout, "%f %f %f\n", m, v, sd); */

    /* sample_sum_array sa; */
    /* sa = init_sample_sum_array(l); */
    /* for (i = 0; i < sa.length; i++) { */
    /*     fprintf(stdout, "%f\n", sa.a[i].sum); */
    /* } */
    /* for (i = 0; i < l; i++) { */
    /*     v1.a[i] = 2; */
    /*     v2.a[i] = 4; */
    /* } */
    /* update_sample_sum_array(&sa, &v1); */
    /* update_sample_sum_array(&sa, &v2); */
    /* for (i = 0; i < sa.length; i++) { */
    /*     fprintf(stdout, "%f\n", sa.a[i].sum_of_squares); */
    /* } */
    /* append_d_array(&v1, 8.0); */
    /* fprintf(stdout, "%d %f\n", v1.length, v1.a[4]); */
    /* c_array c; */
    /* c = init_c_array(4); */
    /* printf("%s %lx %d\n", c.a, strlen(c.a), c.length); */
    /* c.a[0] = 't'; */
    /* c.a[1] = 'e'; */
    /* c.a[2] = 's'; */
    /* c.a[3] = 't'; */
    /* c.a[4] = '\0'; */
    /* printf("%s %lx %d\n", c.a, strlen(c.a), c.length); */
    /* printf("%c\n", c.a[3]); */
    /* if (c.a[4] == '\0') { */
    /*     printf("yes\n"); */
    /* } */
    /* append_c_array(&c, 's'); */
    /* printf("%s %lx %d\n", c.a, strlen(c.a), c.length); */
    /* c_array_2d cc; */
    /* cc = init_c_array_2d(5, 3); */
    /* cc.a[0] = c; */
    /* cc.a[1].a[0] = 'a'; */
    /* cc.a[1].a[1] = 'r'; */
    /* cc.a[1].a[2] = 'e'; */
    /* cc.a[2].a[0] = 'n'; */
    /* cc.a[2].a[1] = 'o'; */
    /* cc.a[2].a[2] = 't'; */
    /* cc.a[2].a[3] = '\0'; */
    /* printf("%s %s %s %d\n", cc.a[0].a, cc.a[1].a, cc.a[2].a, cc.length); */
    /* c_array c2; */
    /* c2 = init_c_array(3); */
    /* c2.a = "fun"; */
    /* append_c_array_2d(&cc, &c2); */
    /* printf("%s %s %s %s %d\n", cc.a[0].a, cc.a[1].a, cc.a[2].a, cc.a[3].a, cc.length); */
    /* s_array strings; */
    /* strings = init_s_array(4); */
    /* for (i = 0; i < cc.length; i++) { */
    /*     strings.a[i] = cc.a[i].a; */
    /* } */
    /* printf("%s %s %s %s %d\n", strings.a[0], strings.a[1], strings.a[2], strings.a[3], strings.length); */
    /* char * x; */
    /* x = "foo"; */
    /* append_s_array(&strings, x); */
    /* printf("%s %s %s %s %s %d\n", strings.a[0], strings.a[1], strings.a[2], strings.a[3], strings.a[4], strings.length); */
    /* printf("%lx\n", sizeof(char)); */
    /* printf("%lx\n", sizeof(int)); */
    /* printf("%lx\n", sizeof(float)); */
    /* printf("%lx\n", sizeof(long)); */
    /* printf("%lx\n", sizeof(double)); */

