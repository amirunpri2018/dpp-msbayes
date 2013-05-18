#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/parsing.c"
#include "test_utils.h"

START_TEST (test_parse_header) {
    /* char * path = "../../test/data/observed_stats_extra_line.txt"; */
    char * path = "data/observed_stats_extra_line.txt";
    c_array * line_buffer;
    s_array * header;
    line_buffer = init_c_array(1023);
    header = init_s_array(1);
    parse_header(path, line_buffer, header);
    ck_assert_msg((header->length == 4), "header length is %d, "
            "expecting %d", header->length, 4);
/* stat.1	stat.2	stat.3	stat.4 */
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("parsing");

    TCase * tc_parse_header = tcase_create("parse_header_test_case");
    tcase_add_test(tc_parse_header, test_parse_header);
    suite_add_tcase(s, tc_parse_header);

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

