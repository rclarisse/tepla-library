#include <assert.h>

#include "rdtsc.h"

#include "../ec_bls509_lcl.h"

#define N 1000

//============================================
//   Feature of Field
//============================================
void test_feature(Field f)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "FieldType: %s\n", field_get_name(f));
    gmp_fprintf(stdout, "   characteristic: %Zx\n", *field_get_char(f));
    gmp_fprintf(stdout, "   order of field: %Zx^%d\n", *field_get_char(f), field_get_degree(f));
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

    char loop[] = "100";

    mpz_t e, exp;

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
    element_set_str(a, "125E40BCE152718F7854FB56EA1D5A7912DDC8A1BD8378AFD3BE21D858E2C652B3CFF5BA6D3262E672458708905887C1BDA1F799795F4150604A58C8C5D65971 EE30691BAEC1B181F0B2FE7D556068A241AE172427584F6AD67670064B7DB948DACD6F54C62BC479749EF3B399F02A122C89817CE7AD8B37A07B389FB514734");
    element_set_str(b, "599EA6C738F27A0B62A150B68C8F9781F8DAEDF5FBE439A7B9DBBF52767A4B0C8A0AF6A94020190A076C4DAD68005D08322C4B838D54951DEA25FCB09761706 9328FA14F9620CE3DC9004C5B4D17E4E62C5B5CA23503F359CCCEAB319687A109C1FF8BDF49B274AD948DC3FDF11FF3C07A0A55585A8FF08627042AF4A7EF79");
    element_set_str(d, "2A2D42955A7CE9431911DAD592566FA66E0AF37C0293DBE1C71DCC9A039DA481A021F65851BC38642F0EF82F51DB9BF5260FF4A43A6F7A19DAFA706168C9DCC 2C03F330B48714A5FE63D7F36E231783EBC748587920A5DD44A34A7B63DD27A350050C1AF93CDCB7513209EC5D54EC1F4DEE565B847D5A35EF1A62737396402");

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
    element_set_str(d, "66CB8F1E4AAB3363677A66F0A6D15F74E22EDD49A005D1FE8C3AAC1FD749C785DC8F0F86A3A879C8AAE829A9019F52BE39151A8746CB30CE526B7173736DED2 1166F4EAF0EE89F6DB65D708BB0FF22FAE3AEC3C62FE241A16884922B156A81A807B325F44D2769453690E963B15D93AAAC3A5DF9A7D81AEBA76336CF63CD04D");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_init_set_str(e, "1B45F16C848B9C476C1D2FF1FD60A0D0C19BBA6F3ECE3CF6C5FCE4FAB7CAD4FF", 16);

    element_pow(c, a, e);
    element_set_str(d, "CEB60811BDE357E879292F045FFA36AB8FFC8F8B8846BD92FCE402A2DD60AB359350A1EFA2BEF0DB446B8FF8F7A68E376C8E66E8B78852DF3C4C4C5712E9953 305A530E53F21B3B4F0673B434D79029406CC6326B8B5EF6809FA16047540451B00760022FBC2118AA886221CBF9E9F9FF5DD4C6BB376F67747FD77669BEC");

    assert(element_cmp(c, d) == 0);

    mpz_clear(e);

    //--------------------
    //  sqr
    //--------------------
    element_sqr(c, a);
    element_mul(d, a, a);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_sqr(c, a);
    }
    t2 = rdtsc();

    printf("element sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);

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

        assert(element_cmp(b, a) == 0);
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

    for (i = 0; i < 100; i++)
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
//   Frobenius Map \phi_p
//============================================
void test_frob(Field f)
{
    int i;
    unsigned long long int t1, t2;
    mpz_t p;
    Element a, b, c;

    mpz_init_set(p, *field_get_char(f));

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
    {
        element_random(a);
        element_pow(b, a, p);

        bls509_fp2_frob_p(c, a);

        assert(element_cmp(b, c) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        bls509_fp2_frob_p(c, a);
    }
    t2 = rdtsc();

    printf("element frob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_clear(p);

    element_clear(a);
    element_clear(b);
    element_clear(c);
}

//============================================
//   i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;
    char a_str[260];

    size_t blen;
    unsigned char b_str[128];

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

    for (i = 0; i < 1000; i++)
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
//  main program
//============================================
int main(void)
{
    Field f;

    // test for methods
    field_init(f, "bls509_fp2");
    test_feature(f);
    test_arithmetic_operation(f);
    test_sqrt(f);
    test_frob(f);
    test_io(f);

    field_clear(f);

    fprintf(stderr, "ok\n");

    return 0;
}
