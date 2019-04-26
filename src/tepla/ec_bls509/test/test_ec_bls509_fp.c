#include <assert.h>

#include "rdtsc.h"

#include "../ec_bls509_lcl.h"

#define N 1000
#define M 100

#define t 80            // security level for SHA
#define MAP_STR "LCIS"  // string for map to point

//============================================
//  Feature of Elliptic Curve
//============================================
void test_feature(const EC_GROUP ec)
{
    fprintf(stdout, "---\n");
    fprintf(stdout, "Elliptic Curve Type: %s\n", curve_get_name(ec));
    fprintf(stdout, "   Y^2 = X^3 + aX + b\n");
    fprintf(stdout, "   a: ");
    element_print(ec->a);
    fprintf(stdout, "   b: ");
    element_print(ec->b);
    gmp_fprintf(stdout, "   field order: %Zx\n", ec->field->order);
    fprintf(stdout, "   generator of curve: ");
    point_print(ec->generator);
    gmp_fprintf(stdout, "   order: %Zx\n", ec->order);
    gmp_fprintf(stdout, "   trace: %Zx\n", ec->trace);
    gmp_fprintf(stdout, "   cofactor: %Zx\n", ec->cofactor);
    fprintf(stdout, "---\n");
}

//============================================
//  test for arithmetic operations of EC
//============================================
void test_arithmetic_operation(const EC_GROUP ec)
{
    int i, j;
    unsigned long long int t1, t2;

    mpz_t scalar;

    EC_POINT x, y, z, w;

    //-------------------
    //  init
    //-------------------
    point_init(x, ec);
    point_init(y, ec);
    point_init(z, ec);
    point_init(w, ec);

    //-------------------
    //  random
    //-------------------
    for (i = 0; i < 100; i++)
    {
        point_random(x);
        point_random(y);

        assert(point_is_on_curve(x));
        assert(point_is_on_curve(y));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_random(x);
    }
    t2 = rdtsc();

    printf("point random: %.2lf [clock]\n", (double)(t2 - t1) / M);

    //-------------------
    //  add/dob
    //-------------------
    point_add(w, w, x);
    point_add(z, y, z);

    assert(point_cmp(w, x) == 0);
    assert(point_cmp(z, y) == 0);

    point_set_infinity(w);
    point_dob(z, w);

    assert(point_is_infinity(z));

    point_set_str(x, "[5DB8A39548ED54E754C8E2CE6A8539BE76A0BDC17C81B16EE9A00C1286BC90470079BB90B9631D2B48B7176642B4E1566947AE4F2E4BDBB4F28F17D84DBB36E,612BE02CD0297191525BBA92F8389CD7F7E99B2CD82FB58961B9D8DBDDFAC06E3A089FB063686C6414F93E77F47BF8080C657E85E4BAA8032F5D3E5FC0EBDAB]");
    point_set_str(y, "[18E88C42C189FBB6C3D9A07634234BB787B9F2F7EBED64442AE918C2E6A01BC42E97584C8FFD41B975447B8C33A753DA880F0FC59DCDADD22411832EA448E97,B4C0C56F18D1FF7C66B8CEC4648837DB8B40C9E6D6140A8AB23EABF477BEB9AF97E79CA96EDD3A094722F63E01246AA70CDAE26FF69322A5F5EFFA813A0B44E]");
    point_set_str(w, "[149B9475C7BA8E86AFBD84C5BAE8CADF68543F99B5C49E1A4FD330C06AFAF0CA7A44B0EBA9A3E480C0F9ABD02D67DF6CC88DD14CA5572AA611A98D26A6C2F6C7,14CF0946A6FF694972EA09786944645B71B83262967E962614B2F17EF8AC122255B5BD8DCCBD5446A49F78FE25AF49D3156ACB43A9ED5461AC3FAAC6B6029F8F]");

    point_add(z, x, y);

    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_add(z, x, y);
    }
    t2 = rdtsc();

    printf("point add: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_set_str(w, "[4314DB1F4AC0BF9ADF2087D7304AC100352A0E403F898591E86D8FAF4F262A94BA9379E1BB80DF4CB9426C4EC662F2F363B5B2C6B9D8A217DD6F74D9C9201,EC1A516A2F2EF7BE5F84DAE7589192DD905D60B4AC41CC50EA471A047C47AB6A4AB774B4598E68AD7564C275866069CEA50C20D5500180EB84F9D13A94C54C7]");

    point_dob(z, x);

    assert(point_cmp(z, w) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_dob(z, x);
    }
    t2 = rdtsc();

    printf("point dob: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //------------------
    // neg/sub
    //------------------
    point_neg(z, x);
    point_add(w, x, z);
    point_sub(z, x, x);

    assert(point_is_infinity(w));
    assert(point_is_infinity(z));

    //-------------------
    //  mul
    //-------------------
    mpz_init(scalar);

    // mpz_init_set_str(scalar, "2370FB049D410FBDC023FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16);
    //
    // ec_bls509_fp_mul(z, scalar, x);
    //
    // ec_bls509_fp_point_endomorphism(z, x);
    //
    // assert(point_cmp(y, z) == 0);
    //
    // point_set_str(x, "[8957DD9D913386BDEE43720A14C1AF42B44C072F78516FBEB12D80DA16A52C9,1520FC7A1090AFE300849C46AA12B287F599C54EC411D498D24E1968B20848A0]");
    point_set_str(y, "[7D18B7E1DDB79EB836817E8AD46D2CA21297C4FB4DEF8380605D7665F51354C03D6E77AE704B9F2A09A6202C6124EC0001E8A52FB9C62C3C50782C20EC49B44,7DFE65CA138F593DCE01689D2F6BD34289C1351A301C956550E8F3ACAC1AF29F21AD53C6C978DD612041EE8B2ABD626A59652F057B7636E4AE49093C4FEDAB8]");

    mpz_set_str(scalar, "1C26BA60159C76F9C565E8716605BFBF2881FF4D82B2B4081F6EB20604CB1830", 16);

    point_mul(z, scalar, x);

    assert(point_cmp(z, y) == 0);

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        point_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul endomorphism: %.2lf [clock]\n", (double)(t2 - t1) / M);

    for (i = 1; i < 100; i++)
    {
        point_random(x);
        point_set_infinity(w);
        mpz_set_ui(scalar, i);

        for (j = 0; j < i; j++) {
            point_add(w, w, x);
        }

        point_mul(z, scalar, x);

        assert(point_cmp(z, w) == 0);
    }

    mpz_set(scalar, ec->order);

    char s[260];

    for (i = 0; i < 100; i++)
    {
        point_random(x);

        point_get_str(s, x);
        fprintf(stderr, "\n");
        fprintf(stderr, s);
        fprintf(stderr, "\n");

        ec_bls509_fp_mul(z, scalar, x);

        point_get_str(s, z);
        fprintf(stderr, "\n");
        fprintf(stderr, s);
        fprintf(stderr, "\n");

        assert(point_is_infinity(z));
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        ec_bls509_fp_mul(z, scalar, x);
    }
    t2 = rdtsc();

    printf("point mul normal: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(scalar);

    //-------------------
    //  clear
    //-------------------
    point_clear(x);
    point_clear(y);
    point_clear(z);
    point_clear(w);
}

// //============================================
// //  MAP to POINT test
// //============================================
// void test_map_to_point(const EC_GROUP ec)
// {
//     int i;
//     unsigned long long int t1, t2;
//     EC_POINT P, Q;
//
//     point_init(P, ec);
//     point_init(Q, ec);
//
//     point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);
//
//     assert(point_is_on_curve(P));
//
//     point_map_to_point(Q, MAP_STR, sizeof(MAP_STR), t);
//
//     assert(point_cmp(Q, P) == 0);
//
//     t1 = rdtsc();
//     for (i = 0; i < N; i++) {
//         point_map_to_point(P, MAP_STR, sizeof(MAP_STR), t);
//     }
//     t2 = rdtsc();
//
//     printf("point map to point in 128 security: %.2lf [clock]\n", (double)(t2 - t1) / N);
//
//     point_clear(P);
//     point_clear(Q);
// }

//============================================
//  i/o test of EC
//============================================
void test_io(const EC_GROUP ec)
{
    int i;
    unsigned long long int t1, t2;

    size_t osize;
    unsigned char os[256];

    char str[512];

    EC_POINT P, Q, R;

    point_init(P, ec);
    point_init(Q, ec);
    point_init(R, ec);

    //---------------------
    //  octet string
    //---------------------
    point_set_infinity(R);

    point_to_oct(os, &osize, R);
    point_from_oct(Q, os, osize);

    assert(point_is_infinity(Q));

    for (i = 0; i < 1000; i++)
    {
        point_random(P);

        point_to_oct(os, &osize, P);
        point_from_oct(Q, os, osize);

        assert(point_cmp(P, Q) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_to_oct(os, &osize, P);
    }
    t2 = rdtsc();

    printf("point to octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_from_oct(Q, os, osize);
    }
    t2 = rdtsc();

    printf("point from octet string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //---------------------
    //  string
    //---------------------
    point_set_infinity(R);

    point_get_str(str, R);
    point_set_str(Q, str);

    assert(point_is_infinity(Q));

    for (i = 0; i < 1000; i++)
    {
        point_get_str(str, P);
        point_set_str(Q, str);

        assert(point_cmp(P, Q) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_get_str(str, P);
    }
    t2 = rdtsc();

    printf("point get string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        point_set_str(Q, str);
    }
    t2 = rdtsc();

    printf("point set string: %.2lf [clock]\n", (double)(t2 - t1) / N);

    point_clear(P);
    point_clear(Q);
    point_clear(R);
}

//============================================
//  main program
//============================================
int main(void)
{
    EC_GROUP ec;

    // test for methods
    curve_init(ec, "ec_bls509_fp");
    test_feature(ec);
    test_arithmetic_operation(ec);
    // test_map_to_point(ec);
    test_io(ec);

    curve_clear(ec);

    fprintf(stderr, "ok\n");

    return 0;
}
