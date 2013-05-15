#include <stdlib.h>
#include <check.h>
#include "../src/array_utils.c"

START_TEST (test_init_d_array) {
    d_array * d;
    d = init_d_array(10);
    ck_assert_int_eq(d->capacity, 10);
    ck_assert_int_eq(d->length, 0);
    free_d_array(d);
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("array_utils");
    TCase * tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_init_d_array);
    suite_add_tcase(s, tc_core);
    
    return s;
}

int main(void) {
    int number_failed;
    Suite * s = array_utils_suite();
    SRunner * sr = srunner_create(s);
    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
