#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/parsing.c"
#include "test_utils.h"

START_TEST (test_parse_header) {
    char * path = "data/observed_stats.txt";
    int ret;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    parse_header(path, line_buffer, header);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
}
END_TEST

START_TEST (test_parse_observed_stats_file) {
    char * path = "data/observed_stats.txt";
    int ret;
    double e = 0.000001;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    d_array * stats;
    d_array * expected_stats;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    stats = init_d_array(1);
    expected_stats = init_d_array(1);
    append_d_array(expected_stats, 0.4);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.4);

    parse_observed_stats_file(path, line_buffer, header, stats);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    ck_assert_msg((stats->length == 4), "stats length is %d, "
            "expecting %d", header->length, 4);
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    ret = d_arrays_equal(stats, expected_stats, e);
    ck_assert_msg((ret != 0), "parsed stats are incorrect");
}
END_TEST

START_TEST (test_parse_observed_stats_file_extra_line) {
    char * path = "data/observed_stats_extra_line.txt";
    int ret;
    double e = 0.000001;
    c_array * line_buffer;
    s_array * header;
    s_array * expected_header;
    d_array * stats;
    d_array * expected_stats;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    expected_header = init_s_array(1);
    append_s_array(expected_header, "stat.1");
    append_s_array(expected_header, "stat.2");
    append_s_array(expected_header, "stat.3");
    append_s_array(expected_header, "stat.4");
    stats = init_d_array(1);
    expected_stats = init_d_array(1);
    append_d_array(expected_stats, 0.1);
    append_d_array(expected_stats, 0.2);
    append_d_array(expected_stats, 0.3);
    append_d_array(expected_stats, 0.4);

    parse_observed_stats_file(path, line_buffer, header, stats);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
    ck_assert_msg((stats->length == 4), "stats length is %d, "
            "expecting %d", header->length, 4);
    ret = s_arrays_equal(header, expected_header);
    ck_assert_msg((ret != 0), "parsed header is incorrect");
    ret = d_arrays_equal(stats, expected_stats, e);
    ck_assert_msg((ret != 0), "parsed stats are incorrect");
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("parsing");

    TCase * tc_parse_header = tcase_create("parse_header_test_case");
    tcase_add_test(tc_parse_header, test_parse_header);
    suite_add_tcase(s, tc_parse_header);

    TCase * tc_parse_observed_stats_file = tcase_create(
            "parse_observed_stats_file_test_case");
    tcase_add_test(tc_parse_observed_stats_file,
            test_parse_observed_stats_file);
    tcase_add_test(tc_parse_observed_stats_file,
            test_parse_observed_stats_file_extra_line);
    suite_add_tcase(s, tc_parse_observed_stats_file);

    return s;
}

int main(void) {
    int number_failed;
    Suite * s = array_utils_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_VERBOSE);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

