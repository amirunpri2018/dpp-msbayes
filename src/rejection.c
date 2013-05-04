/**
 * \file rejection.c
 * \brief A collection of functions for Euclidean distance rejection.
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h> // for getopt
#include <string.h>

#define VERSION "0.1"

typedef struct d_array_ {
    double * a;
    int length;
} d_array;

typedef struct c_array_ {
    char * a;
    int length;
} c_array;

typedef struct c_array_2d_ {
    c_array * a;
    int length;
} c_array_2d;

typedef struct s_array_ {
    char ** a;
    int length;
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
} config;

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

void append_d_array(d_array * v, double x) {
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

d_array get_mean_array(const sample_sum_array * s) {
    d_array v;
    v = init_d_array((*s).length);
    int i;
    for (i = 0; i < (*s).length; i++) {
        v.a[i] = get_mean(&(*s).a[i]);
    }
    return v;
}

d_array get_sample_variance_array(const sample_sum_array * s) {
    d_array v;
    v = init_d_array((*s).length);
    int i;
    for (i = 0; i < (*s).length; i++) {
        v.a[i] = get_sample_variance(&(*s).a[i]);
    }
    return v;
}
        
d_array get_std_dev_array(const sample_sum_array * s) {
    d_array v;
    v = init_d_array((*s).length);
    int i;
    for (i = 0; i < (*s).length; i++) {
        v.a[i] = get_std_dev(&(*s).a[i]);
    }
    return v;
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
    fprintf(stderr, "Path(s) to files with simulated draws:\n");
    for (i = 0; i < (*c).sim_paths.length; i++) {
        fprintf(stderr, "\t%s\n", (*c).sim_paths.a[i]);
    }
}

void parse_args(config * conf, int argc, char **argv) {
    int i, j, k;
    char * p;
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

int main(int argc, char **argv){
    config conf;
    conf = init_config();
    parse_args(&conf, argc, argv);

    int i, l = 4;
    d_array v1, v2, means, std_devs;
    v1 = init_d_array(l);
    v2 = init_d_array(l);
    means = init_d_array(l);
    std_devs = init_d_array(l);
    for (i = 0; i < l; i++) {
        v1.a[i] = 2;
        v2.a[i] = 4;
        means.a[i] = 3;
        std_devs.a[i] = 2;
    }
    double d1 = get_euclidean_distance(&v1, &v2);
    standardize_vector(&v1, &means, &std_devs);
    standardize_vector(&v2, &means, &std_devs);
    double d2 = get_euclidean_distance(&v1, &v2);
    fprintf(stdout, "%f\n%f\n", d1, d2);
    /* d_array v3 = init_d_array(5); */
    /* double d3 = get_euclidean_distance(v1, v3); */

    sample_sum s;
    s = init_sample_sum();
    update_sample_sum(&s, 2.0);
    update_sample_sum(&s, 4.0);
    update_sample_sum(&s, 6.0);
    double m, v, sd;
    m = get_mean(&s);
    v = get_sample_variance(&s);
    sd = get_std_dev(&s);
    fprintf(stdout, "%f %f %f\n", m, v, sd);

    sample_sum_array sa;
    sa = init_sample_sum_array(l);
    for (i = 0; i < sa.length; i++) {
        fprintf(stdout, "%f\n", sa.a[i].sum);
    }
    for (i = 0; i < l; i++) {
        v1.a[i] = 2;
        v2.a[i] = 4;
    }
    update_sample_sum_array(&sa, &v1);
    update_sample_sum_array(&sa, &v2);
    for (i = 0; i < sa.length; i++) {
        fprintf(stdout, "%f\n", sa.a[i].sum_of_squares);
    }
    append_d_array(&v1, 8.0);
    fprintf(stdout, "%d %f\n", v1.length, v1.a[4]);
    c_array c;
    c = init_c_array(4);
    printf("%s %lx %d\n", c.a, strlen(c.a), c.length);
    c.a[0] = 't';
    c.a[1] = 'e';
    c.a[2] = 's';
    c.a[3] = 't';
    c.a[4] = '\0';
    printf("%s %lx %d\n", c.a, strlen(c.a), c.length);
    printf("%c\n", c.a[3]);
    if (c.a[4] == '\0') {
        printf("yes\n");
    }
    append_c_array(&c, 's');
    printf("%s %lx %d\n", c.a, strlen(c.a), c.length);
    c_array_2d cc;
    cc = init_c_array_2d(5, 3);
    cc.a[0] = c;
    cc.a[1].a[0] = 'a';
    cc.a[1].a[1] = 'r';
    cc.a[1].a[2] = 'e';
    cc.a[2].a[0] = 'n';
    cc.a[2].a[1] = 'o';
    cc.a[2].a[2] = 't';
    cc.a[2].a[3] = '\0';
    printf("%s %s %s %d\n", cc.a[0].a, cc.a[1].a, cc.a[2].a, cc.length);
    c_array c2;
    c2 = init_c_array(3);
    c2.a = "fun";
    append_c_array_2d(&cc, &c2);
    printf("%s %s %s %s %d\n", cc.a[0].a, cc.a[1].a, cc.a[2].a, cc.a[3].a, cc.length);
    s_array strings;
    strings = init_s_array(4);
    for (i = 0; i < cc.length; i++) {
        strings.a[i] = cc.a[i].a;
    }
    printf("%s %s %s %s %d\n", strings.a[0], strings.a[1], strings.a[2], strings.a[3], strings.length);
    char * x;
    x = "foo";
    append_s_array(&strings, x);
    printf("%s %s %s %s %s %d\n", strings.a[0], strings.a[1], strings.a[2], strings.a[3], strings.a[4], strings.length);
    printf("%lx\n", sizeof(char));
    printf("%lx\n", sizeof(int));
    printf("%lx\n", sizeof(float));
    printf("%lx\n", sizeof(long));
    printf("%lx\n", sizeof(double));
}

