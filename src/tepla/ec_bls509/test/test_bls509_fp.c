#include <assert.h>

#include "rdtsc.h"

#include "../ec_bls509_lcl.h"

#define N 10000

//============================================
//   Feature of Field
//============================================
void test_feature(Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "FieldType: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   order of field: %Zx\n", f->order);
    fprintf(stdout, "---\n");
}

//============================================
//   test for arithmetic operations
//============================================
void test_arithmetic_operation(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    char loop[] = "1000";

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "11E220FD094D03865735F8C9A820DA856D497BC8E3510C043B48437B188B25B86C6FAF49CB978D5005609BEEBCCDF0168C85168DF20C2039B8EA4B74D5D1A1DA");
    element_set_str(b, "A7CFDAC68D641AC212FC958E76BDD30F942E6FE1C405CCA175285DB91438553274B21A62859DD80B4CC3F13A4AA7CDFB990AE19EB8A4895368FB4829C8289C0");
    element_set_str(d, "709C7A972E97A967B77CF6D95CBCABF9B019A7DA278EA421FB0C852C9BE1A50314C4B3077D8C9DFEA617EA1EFBD992357B207A06F08D5CE4E3CEE69B99458EF");

    element_add(c, a, b);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_add(c, a, b);
    }
    t2 = rdtsc();

    printf("element add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  sub
    //--------------------
    element_set(d, c);
    element_sub(c, c, d);

    assert(element_is_zero(c));

    //--------------------
    //  mul
    //--------------------
    element_mul(c, a, b);
    element_set_str(d, "4947E5C2EAB554944887C0BE2DC937631AD3199F23B56664A517D5BB7A8F33B53FFC122E07278B3A152E64DF9AC8F27626A953456CEA4A132D1D4FF236141DB");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  random
    //--------------------
    element_random(a);
    element_random(b);

    //--------------------
    //  inv
    //--------------------
    element_mul(c, a, b);
    element_inv(b, b);
    element_mul(c, c, b);
    element_inv(d, a);
    element_mul(d, a, d);

    assert(element_cmp(c, a) == 0);
    assert(element_is_one(d));

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_inv(b, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------

    mpz_init_set_str(exp, loop, 10);

    element_set_one(b);

    for (i = 0; i < atoi(loop); i++) {
        element_mul(b, b, a);
    }

    element_pow(c, a, exp);

    assert(element_cmp(b, c) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, exp);

        assert(element_cmp(a, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//   test for sqrt
//============================================
void test_sqrt(Field f)
{
    int i;
    unsigned long long int t1, t2;
    Element a, b, c, d;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    for (i = 0; i < 1000; i++)
    {
        element_random(a);
        element_sqr(b, a);

        assert(element_is_sqr(b));

        element_sqrt(c, b);
        element_sqr(d, c);

        assert(element_cmp(d, b) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_is_sqr(b);
    }
    t2 = rdtsc();


    printf("element is sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqrt(c, b);
    }
    t2 = rdtsc();

    printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//  i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;

    char a_str[256];

    size_t blen;
    unsigned char b_str[256];

    Element a, b, c;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 1000; i++)
    {
        element_random(a);

        element_get_str(a_str, a);
        element_set_str(c, a_str);

        assert(element_cmp(a, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_get_str(a_str, a);
    }
    t2 = rdtsc();

    printf("element get string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_set_str(c, a_str);
    }
    t2 = rdtsc();

    printf("element set string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    for (i = 0; i < 10000; i++)
    {
        element_random(b);

        element_to_oct(b_str, &blen, b);
        element_from_oct(c, b_str, blen);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_to_oct(b_str, &blen, b);
    }
    t2 = rdtsc();

    printf("element to octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_from_oct(c, b_str, blen);
    }
    t2 = rdtsc();

    printf("element from octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
// main program
//============================================
int main(void)
{
    Field f;

    field_init(f, "bls509_fp");
    test_feature(f);
    test_arithmetic_operation(f);
    test_sqrt(f);
    test_io(f);

    field_clear(f);

    fprintf(stderr, "ok\n");

    return 0;
}
