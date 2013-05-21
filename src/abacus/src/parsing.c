/**
 * @file        parsing.c
 * @authors     Jamie Oaks
 * @package     ABACUS (Approximate BAyesian C UtilitieS)
 * @brief       A collection of parsing functions.
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

#include "parsing.h"

void parse_header(const char * path, c_array * line_buffer, s_array * header) {
    FILE * f;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no lines in %s", path);
        exit(1);
    }
    split_str((*line_buffer).a, header, 0);
    fclose(f);
}

void parse_observed_stats_file(const char * path, c_array * line_buffer,
        s_array * header, d_array * stats) {
    FILE * f;
    int column_idx, n;
    char * ptr;
    double fmatch;
    (*header).length = 0;
    (*stats).length = 0;
    if ((f = fopen(path, "r")) == NULL) {
        perror(path);
        exit(1);
    }
    // parse header
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no header in %s\n", path);
        exit(1);
    }
    split_str((*line_buffer).a, header, 0);
    // parse stats
    if ((fgets((*line_buffer).a, (((*line_buffer).capacity) - 1), f)) == NULL) {
        fprintf(stderr, "ERROR: found no stats in %s\n", path);
        exit(1);
    }
    column_idx = 0;
    ptr = (*line_buffer).a;
    while(*ptr) {
        if ((sscanf(ptr, "%lf%n", &fmatch, &n)) == 1) {
            ptr += n;
            append_d_array(stats, fmatch);
            column_idx++;
        }
        ++ptr;
    }
    fclose(f);
    if ((*header).length != (*stats).length) {
        fprintf(stderr, "ERROR: found %d column headers, but %d stats in "
                "file %s\n", (*header).length, (*stats).length, path);
        exit(1);
    }
}
