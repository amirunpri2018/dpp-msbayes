#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/array_utils.h"
#include "../src/array_utils.c"
#include "../src/stats_utils.h"
#include "../src/stats_utils.c"

START_TEST (test_init_free_sample_sum) {
    sample_sum * ss;
    ss = init_sample_sum();
    ck_assert_msg((ss->n == 0), "`sample_sum.n` did not initiate to 0");
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("stats_utils");

    TCase * tc_sample_sum = tcase_create("sample_sum_test_case");
    tcase_add_test(tc_sample_sum, test_init_free_sample_sum);
    suite_add_tcase(s, tc_sample_sum);

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

