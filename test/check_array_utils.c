#include <stdlib.h>
#include <check.h>
#include "../src/array_utils.c"

START_TEST (test_init_free_d_array) {
    int i, size;
    d_array * v;
    size = 10;
    v = init_d_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`d_array.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`d_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `d_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(double)), "Pointer to array"
            "of `d_array` is not pointing to double");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i] == 0), "Element %d of %d was not initialized "
            "to zero", i, size);
    }
    free_d_array(v);
}
END_TEST

START_TEST (test_append_d_array) {
    int i, size;
    d_array * v;
    size = 10;
    v = init_d_array(size);
    ck_assert_int_eq(v->length, 0);
    append_d_array(v, 0.1);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ck_assert_msg((v->a[0] == 0.1), "`d_array` element is %lf, expecting %lf",
            v->a[0], 0.1);
    append_d_array(v, 1.9);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 2);
    ck_assert_msg((v->a[1] == 1.9), "`d_array` element is %lf, expecting %lf",
            v->a[0], 1.9);
}
END_TEST


START_TEST (test_init_free_i_array) {
    int i, size;
    i_array * v;
    size = 10;
    v = init_i_array(size);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 0);
    ck_assert_msg((sizeof(v->length) == sizeof(int)), "`i_array.length` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->capacity) == sizeof(int)), "`i_array.capacity` "
            "attribute is not an int");
    ck_assert_msg((sizeof(v->a) == sizeof(long)), "Array attribute "
            "of `i_array` is not a pointer");
    ck_assert_msg((sizeof(*v->a) == sizeof(int)), "Pointer to array"
            "of `i_array` is not pointing to int");
    for (i = 0; i < v->capacity; i++) {
        ck_assert_msg((v->a[i] == 0), "Element %d of %d was not initialized "
            "to zero", i, size);
    }
    free_i_array(v);
}
END_TEST

START_TEST (test_append_i_array) {
    int i, size;
    i_array * v;
    size = 10;
    v = init_i_array(size);
    ck_assert_int_eq(v->length, 0);
    append_i_array(v, 1);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 1);
    ck_assert_msg((v->a[0] == 1), "`i_array` element is %lf, expecting %lf",
            v->a[0], 1);
    append_i_array(v, -99);
    ck_assert_int_eq(v->capacity, size);
    ck_assert_int_eq(v->length, 2);
    ck_assert_msg((v->a[1] == -99), "`i_array` element is %lf, expecting %lf",
            v->a[0], -99);
}
END_TEST

Suite * array_utils_suite(void) {
    Suite * s = suite_create("array_utils");

    TCase * tc_d_array = tcase_create("d_array_test_case");
    tcase_add_test(tc_d_array, test_init_free_d_array);
    tcase_add_test(tc_d_array, test_append_d_array);
    suite_add_tcase(s, tc_d_array);

    TCase * tc_i_array = tcase_create("i_array_test_case");
    tcase_add_test(tc_i_array, test_init_free_i_array);
    tcase_add_test(tc_i_array, test_append_i_array);
    suite_add_tcase(s, tc_i_array);
    
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

