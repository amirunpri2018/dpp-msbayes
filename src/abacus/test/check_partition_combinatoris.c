#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/partition_combinatorics.c"
#include "test_utils.h"


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

