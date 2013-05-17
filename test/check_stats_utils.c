#include <stdlib.h>
#include <check.h>
#include <signal.h>
#include "../src/stats_utils.c"
#include "test_utils.h"

START_TEST (test_init_free_sample_sum) {
    sample_sum * ss;
    ss = init_sample_sum();
    ck_assert_msg((ss->n == 0), "`sample_sum.n` did not initiate to 0");
    ck_assert_msg((ss->sum == 0.0), "`sample_sum.sum` did not initiate to 0");
    ck_assert_msg((ss->sum_of_squares == 0.0), "`sample_sum.sum_of_squares` "
            "did not initiate to 0");
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_update_sample_sum) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 1.0);
    ck_assert_msg((ss->n == 1), "`sample_sum.n` is %d, expecting %d", ss->n, 1);
    ck_assert_msg(almost_equal(ss->sum, 1.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 1.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 1.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 1.0);

    update_sample_sum(ss, 2.0);
    ck_assert_msg((ss->n == 2), "`sample_sum.n` is %d, expecting %d", ss->n, 2);
    ck_assert_msg(almost_equal(ss->sum, 3.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 3.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 5.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 5.0);

    update_sample_sum(ss, 3.0);
    ck_assert_msg((ss->n == 3), "`sample_sum.n` is %d, expecting %d", ss->n, 2);
    ck_assert_msg(almost_equal(ss->sum, 6.0, e), "`sample_sum.sum` is %lf, "
            "expecting %lf", ss->sum, 6.0);
    ck_assert_msg(almost_equal(ss->sum_of_squares, 14.0, e),
            "`sample_sum.sum_of_squares is %lf, expecting %lf", ss->sum, 14.0);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_get_mean_n0) {
    sample_sum * ss;
    ss = init_sample_sum();
    get_mean(ss);
}
END_TEST

START_TEST (test_get_mean) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    ck_assert_msg(almost_equal(get_mean(ss), 100.0, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 100.0);
    update_sample_sum(ss, 101.0);
    ck_assert_msg(almost_equal(get_mean(ss), 100.5, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 100.5);
    update_sample_sum(ss, 117.0);
    ck_assert_msg(almost_equal(get_mean(ss), 106.0, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), 106.0);
    update_sample_sum(ss, -500.0);
    ck_assert_msg(almost_equal(get_mean(ss), -45.5, e), "mean is %lf, "
            "expecting %lf", get_mean(ss), -45.5);
    free_sample_sum(ss);
}
END_TEST

START_TEST (test_get_sample_variance_n0) {
    sample_sum * ss;
    ss = init_sample_sum();
    get_sample_variance(ss);
}
END_TEST

START_TEST (test_get_sample_variance_n1) {
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    get_sample_variance(ss);
}
END_TEST

START_TEST (test_get_sample_variance) {
    double e = 0.000001;
    sample_sum * ss;
    ss = init_sample_sum();
    update_sample_sum(ss, 100.0);
    update_sample_sum(ss, 101.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 0.5, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss), 0.5);
    update_sample_sum(ss, 117.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 91.0, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss), 91.0);
    update_sample_sum(ss, -500.0);
    ck_assert_msg(almost_equal(get_sample_variance(ss), 91869.66666667, e),
            "variance is %lf, expecting %lf", get_sample_variance(ss),
            91869.666666667);
    free_sample_sum(ss);
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("stats_utils");

    TCase * tc_sample_sum = tcase_create("sample_sum_test_case");
    tcase_add_test(tc_sample_sum, test_init_free_sample_sum);
    tcase_add_test(tc_sample_sum, test_update_sample_sum);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_mean_n0, SIGABRT);
    tcase_add_test(tc_sample_sum, test_get_mean);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_sample_variance_n0,
            SIGABRT);
    tcase_add_test_raise_signal(tc_sample_sum, test_get_sample_variance_n1,
            SIGABRT);
    tcase_add_test(tc_sample_sum, test_get_sample_variance);
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

