//----------------------------------------------------
//  Header file for ec_bls509
//----------------------------------------------------
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <tepla/ec.h>
#include <tepla/hash.h>

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

//---------------------------------------------------
// Finite Field (bls509) ID
//---------------------------------------------------
typedef enum
{
    bls509_fp,
    bls509_fp2,
    bls509_fp4,
    bls509_fp12,
    bls509_fp24,

} BLS509_FieldType;

//---------------------------------------------------
//  precomputation values for frobenius map
//---------------------------------------------------
typedef struct ec_field_precomp_frob_st
{
    size_t glen1;
    size_t glen2;
    size_t glen3;

    Element *gamma1;
    Element *gamma2;
    Element *gamma3;

} *field_precomp_frob_p;

//---------------------------------------------------
// structure for precomputation values
//---------------------------------------------------
typedef struct ec_field_precomp_st
{
    field_precomp_frob_p  pf;

} *field_precomp_p;

//---------------------------------------------------
// Elliptic Curve (BLS509) ID
//---------------------------------------------------
typedef enum
{
    ec_bls509_fp,
    ec_bls509_fp4,

} BLS509_CurveType;

//---------------------------------------------------
// structure for Elliptic Curve values
//---------------------------------------------------
typedef struct ec_bls509_fp_ec_data_st
{
    Element beta;

    mpz_t n;
    mpz_t n2;
    mpz_t a1, a2;
    mpz_t b1, b2;

} *ec_data_fp;

//---------------------------------------------------
// structure for Elliptic Curve values
//---------------------------------------------------
typedef struct ec_bls509_fp4_ec_data_st
{
    mpz_t _6x;
    mpz_t _6x2;

    Element vfrobx, vfroby;
    Element vfrobx2, vfroby2;
    Element vfrobx3, vfroby3;

} *ec_data_fp4;

//---------------------------------------------------
//  structure for precomputation value
//---------------------------------------------------
typedef struct ec_pairing_precomp_st
{
    size_t slen;
    int *si;

    size_t tlen; // for calculating f^t
    int *ti;     // for calculating f^t

} *pairing_precomp_p;

//----------------------------------------------
// declaration function of field bls509_fp
//----------------------------------------------
void bls509_fp_init(Element x);
void bls509_fp_clear(Element x);
void bls509_fp_set(Element x, const Element y);
void bls509_fp_set_str(Element x, const char *str);
void bls509_fp_get_str(char *str, const Element x);
void bls509_fp_set_zero(Element x);
void bls509_fp_set_one(Element x);
void bls509_fp_mod(Element z, const Element x);
void bls509_fp_add(Element z, const Element x, const Element y);
void bls509_fp_addn(Element z, const Element x, const Element y);
void bls509_fp_addp(Element z, const Element x);
void bls509_fp_add_one(Element z, const Element x);
void bls509_fp_dob(Element z, const Element x);
void bls509_fp_tri(Element z, const Element x);
void bls509_fp_neg(Element z, const Element x);
void bls509_fp_sub(Element z, const Element x, const Element y);
void bls509_fp_subn(Element z, const Element x, const Element y);
void bls509_fp_mul(Element z, const Element x, const Element y);
void bls509_fp_muln(Element z, const Element x, const Element y);
void bls509_fp_mulc(Element z, const Element x, const mpz_t c);
void bls509_fp_div2(Element z, const Element x);
void bls509_fp_sqr(Element z, const Element x);
void bls509_fp_inv(Element z, const Element x);
void bls509_fp_pow(Element z, const Element x, const mpz_t exp);
int  bls509_fp_sqrt(Element z, const Element x);
void bls509_fp_OP1_1(Element z, const Element x);
void bls509_fp_OP1_2(Element z, const Element x);
void bls509_fp_OP2(Element z, const Element x);
int  bls509_fp_is_zero(const Element x);
int  bls509_fp_is_one(const Element x);
int  bls509_fp_is_sqr(const Element x);
int  bls509_fp_is_sqr_general(const Element x);
int  bls509_fp_cmp(const Element x, const Element y);
void bls509_fp_random(Element z);
void bls509_fp_to_oct(unsigned char *os, size_t *size, const Element x);
void bls509_fp_from_oct(Element z, const unsigned char *os, const size_t size);
void bls509_fp_FE2IP(mpz_t dst, const Element x);

//----------------------------------------------
// declaration function of field bls509_fp2
//----------------------------------------------
void bls509_fp2_init(Element x);
void bls509_fp2_clear(Element x);
void bls509_fp2_set(Element x, const Element y);
void bls509_fp2_set_fp(Element z, const Element x, const Element y);
void bls509_fp2_set_str(Element x, const char *s);
void bls509_fp2_get_str(char *s, const Element x);
void bls509_fp2_set_zero(Element x);
void bls509_fp2_set_one(Element x);
void bls509_fp2_add(Element z, const Element x, const Element y);
void bls509_fp2_add_one(Element z, const Element x);
void bls509_fp2_addn(Element z, const Element x, const Element y);
void bls509_fp2_dob(Element z, const Element x);
void bls509_fp2_tri(Element z, const Element x);
void bls509_fp2_neg(Element z, const Element x);
void bls509_fp2_sub(Element z, const Element x, const Element y);
void bls509_fp2_subn(Element z, const Element x, const Element y);
void bls509_fp2_mul(Element z, const Element x, const Element y);
void bls509_fp2_muln(Element z, const Element x, const Element y);
void bls509_fp2_mul_p(Element z, const Element x, const Element y);
void bls509_fp2_mul_c(Element z, const Element x, const mpz_t c);
void bls509_fp2_div_2(Element z, const Element x);
void bls509_fp2_sqr(Element z, const Element x);
void bls509_fp2_sqrn(Element z, const Element x);
void bls509_fp2_xi_mul(Element z, const Element x);
// void bls509_fp2_xi_mul_inv(Element z, const Element x);
void bls509_fp2_inv(Element z, const Element x);
void bls509_fp2_pow(Element z, const Element x, const mpz_t exp);
int  bls509_fp2_sqrt(Element z, const Element x);
void bls509_fp2_mod(Element z, const Element x);
void bls509_fp2_OP1_1(Element z, const Element x);
void bls509_fp2_OP1_2(Element z, const Element x);
void bls509_fp2_OP2(Element z, const Element x);
void bls509_fp2_frob_p(Element z, const Element x);
void bls509_fp2_conj(Element z, const Element x);
int  bls509_fp2_is_zero(const Element x);
int  bls509_fp2_is_one(const Element x);
int  bls509_fp2_is_sqr(const Element x);
int  bls509_fp2_cmp(const Element x, const Element y);
void bls509_fp2_precomp(Field f);
void bls509_fp2_random(Element z);
void bls509_fp2_to_oct(unsigned char *os, size_t *size, const Element x);
void bls509_fp2_from_oct(Element z, const unsigned char *os, const size_t size);

// //----------------------------------------------
// // declaration function of field bls509_fp4
// //----------------------------------------------
void bls509_fp4_init(Element x);
void bls509_fp4_clear(Element x);
void bls509_fp4_set(Element x, const Element y);
void bls509_fp4_set_fp2(Element z, const Element x, const Element y);
void bls509_fp4_set_str(Element x, const char *s);
void bls509_fp4_get_str(char *s, const Element x);
void bls509_fp4_set_zero(Element x);
void bls509_fp4_set_one(Element x);
void bls509_fp4_add(Element z, const Element x, const Element y);
void bls509_fp4_addn(Element z, const Element x, const Element y);
void bls509_fp4_dob(Element z, const Element x);
void bls509_fp4_tri(Element z, const Element x);
void bls509_fp4_neg(Element z, const Element x);
void bls509_fp4_sub(Element z, const Element x, const Element y);
void bls509_fp4_subn(Element z, const Element x, const Element y);
void bls509_fp4_mul(Element z, const Element x, const Element y);
void bls509_fp4_muln(Element z, const Element x, const Element y);
void bls509_fp4_beta_mul(Element z, const Element x);

void bls509_fp4_mul_fp2(Element z, const Element x, const Element y);
void bls509_fp4_mul_fp2_2(Element z, const Element y, const Element x1, const Element x2);

void bls509_fp4_mul_fp2_3(Element z, const Element x, const Element y);
void bls509_fp4_mul_fp2_4(Element z, const Element y, const Element x1, const Element x2);

void bls509_fp4_pow(Element z, const Element x, const mpz_t exp);
void bls509_fp4_conj(Element z, const Element x);
void bls509_fp4_sqr(Element z, const Element x);
int  bls509_fp4_sqrt(Element z, const Element x);
void bls509_fp4_inv(Element z, const Element x);
void bls509_fp4_mod(Element z, const Element x);
void bls509_fp4_OP1_1(Element z, const Element x);
void bls509_fp4_OP1_2(Element z, const Element x);
void bls509_fp4_OP2(Element z, const Element x);
void bls509_fp4_frob_p(Element z, const Element x);
int  bls509_fp4_is_zero(const Element x);
int  bls509_fp4_is_one(const Element x);
int  bls509_fp4_is_sqr(const Element x);
int  bls509_fp4_cmp(const Element x, const Element y);
void bls509_fp4_precomp(Field f);
void bls509_fp4_precomp_for_pairing_init(Field f);
void bls509_fp4_random(Element z);
void bls509_fp4_to_oct(unsigned char *os, size_t *size, const Element x);
void bls509_fp4_from_oct(Element z, const unsigned char *os, const size_t size);

// //----------------------------------------------
// // declaration function of field bls509_fp12
// //----------------------------------------------
// void bls509_fp12_init(Element x);
// void bls509_fp12_clear(Element x);
// void bls509_fp12_set(Element x, const Element y);
// void bls509_fp12_set_fp6(Element z, const Element x, const Element y);
// void bls509_fp12_set_str(Element x, const char *s);
// void bls509_fp12_get_str(char *s, const Element x);
// void bls509_fp12_set_zero(Element x);
// void bls509_fp12_set_one(Element x);
// void bls509_fp12_add(Element z, const Element x, const Element y);
// void bls509_fp12_dob(Element z, const Element x);
// void bls509_fp12_tri(Element z, const Element x);
// void bls509_fp12_neg(Element z, const Element x);
// void bls509_fp12_sub(Element z, const Element x, const Element y);
// void bls509_fp12_mul(Element z, const Element x, const Element y);
// void bls509_fp12_mul_L(Element z, Element x0, Element x1, Element x2);
// void bls509_fp12_mul_L2(Element z, Element x0, Element x1, Element x2);
// void bls509_fp12_sqr(Element z, const Element x);
// int  bls509_fp12_sqrt(Element z, const Element x);
// void bls509_fp12_inv(Element z, const Element x);
// void bls509_fp12_pow(Element z, const Element x, const mpz_t exp);
// void bls509_fp12_pow_naf(Element z, const Element x, const mpz_t exp);
// void bls509_fp12_frob_p(Element z, const Element x);
// void bls509_fp12_frob_p2(Element z, const Element x);
// void bls509_fp12_frob_p3(Element z, const Element x);
// void bls509_fp12_conj(Element z, const Element x);
// void bls509_fp12_pow_forpairing(Element z, const Element x, const int *t, int tlen);
// void bls509_fp12_sqr_forpairing_karabina(Element z, const Element x);
// void bls509_fp12_pow_forpairing_karabina(Element z, const Element x, const int *t, int tlen);
// void bls509_fp12_decompose_forpairing_karabina(Element z, const Element x);
// void bls509_fp12_sqr_forpairing_beuchat(Element z, const Element x);
// void bls509_fp12_pow_forpairing_beuchat(Element z, const Element x, const int *t, int tlen);
//
// int  bls509_fp12_is_zero(const Element x);
// int  bls509_fp12_is_one(const Element x);
// int  bls509_fp12_is_sqr(const Element x);
// int  bls509_fp12_cmp(const Element x, const Element y);
// void bls509_fp12_precomp(Field f);
// void bls509_fp12_precomp_for_pairing_init(Field f);
// void bls509_fp12_random(Element z);
// void bls509_fp12_to_oct(unsigned char *os, size_t *size, const Element x);
// void bls509_fp12_from_oct(Element z, const unsigned char *os, const size_t size);
//
// //----------------------------------------------
// // declaration function of field bls509_fp24
// //----------------------------------------------
// void bls509_fp24_init(Element x);
// void bls509_fp24_clear(Element x);
// void bls509_fp24_set(Element x, const Element y);
// void bls509_fp24_set_fp6(Element z, const Element x, const Element y);
// void bls509_fp24_set_str(Element x, const char *s);
// void bls509_fp24_get_str(char *s, const Element x);
// void bls509_fp24_set_zero(Element x);
// void bls509_fp24_set_one(Element x);
// void bls509_fp24_add(Element z, const Element x, const Element y);
// void bls509_fp24_dob(Element z, const Element x);
// void bls509_fp24_tri(Element z, const Element x);
// void bls509_fp24_neg(Element z, const Element x);
// void bls509_fp24_sub(Element z, const Element x, const Element y);
// void bls509_fp24_mul(Element z, const Element x, const Element y);
// void bls509_fp24_mul_L(Element z, Element x0, Element x1, Element x2);
// void bls509_fp24_mul_L2(Element z, Element x0, Element x1, Element x2);
// void bls509_fp24_sqr(Element z, const Element x);
// int  bls509_fp24_sqrt(Element z, const Element x);
// void bls509_fp24_inv(Element z, const Element x);
// void bls509_fp24_pow(Element z, const Element x, const mpz_t exp);
// void bls509_fp24_pow_naf(Element z, const Element x, const mpz_t exp);
// void bls509_fp24_frob_p(Element z, const Element x);
// void bls509_fp24_frob_p2(Element z, const Element x);
// void bls509_fp24_frob_p3(Element z, const Element x);
// void bls509_fp24_conj(Element z, const Element x);
// void bls509_fp24_pow_forpairing(Element z, const Element x, const int *t, int tlen);
// void bls509_fp24_sqr_forpairing_karabina(Element z, const Element x);
// void bls509_fp24_pow_forpairing_karabina(Element z, const Element x, const int *t, int tlen);
// void bls509_fp24_decompose_forpairing_karabina(Element z, const Element x);
// void bls509_fp24_sqr_forpairing_beuchat(Element z, const Element x);
// void bls509_fp24_pow_forpairing_beuchat(Element z, const Element x, const int *t, int tlen);
//
// int  bls509_fp24_is_zero(const Element x);
// int  bls509_fp24_is_one(const Element x);
// int  bls509_fp24_is_sqr(const Element x);
// int  bls509_fp24_cmp(const Element x, const Element y);
// void bls509_fp24_precomp(Field f);
// void bls509_fp24_precomp_for_pairing_init(Field f);
// void bls509_fp24_random(Element z);
// void bls509_fp24_to_oct(unsigned char *os, size_t *size, const Element x);
// void bls509_fp24_from_oct(Element z, const unsigned char *os, const size_t size);


//----------------------------------------------
// declaration function of elliptic curve
//----------------------------------------------
void ec_bls509_fp_point_init(EC_POINT p);
void ec_bls509_fp_point_clear(EC_POINT p);
void ec_bls509_fp_point_set(EC_POINT z, const EC_POINT x);
void ec_bls509_fp_point_set_str(EC_POINT z, const char* s);
void ec_bls509_fp_point_set_infinity(EC_POINT z);
void ec_bls509_fp_point_get_str(char *s, const EC_POINT z);
void ec_bls509_fp_add(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bls509_fp_dob(EC_POINT z, const EC_POINT x);
void ec_bls509_fp_add_formul(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bls509_fp_dob_formul(EC_POINT z, const EC_POINT x);
void ec_bls509_fp_neg(EC_POINT z, const EC_POINT x);
void ec_bls509_fp_sub(EC_POINT z, const EC_POINT x, const EC_POINT y);
void ec_bls509_fp_mul_affine(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bls509_fp_mul(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bls509_fp_mul_naf(EC_POINT z, const mpz_t s, const EC_POINT x);
void ec_bls509_fp_mul_end(EC_POINT z, const mpz_t s, const EC_POINT x);
int  ec_bls509_fp_is_infinity(const EC_POINT P);
int  ec_bls509_fp_is_on_curve(const EC_POINT P);
int  ec_bls509_fp_cmp(const EC_POINT x, const EC_POINT y);
void ec_bls509_fp_make_affine(EC_POINT z, const EC_POINT x);
void ec_bls509_fp_map_to_point(EC_POINT z, const char *s, size_t slen, int t);
void ec_bls509_fp_point_endomorphism(EC_POINT Q, const EC_POINT P);
void ec_bls509_fp_random(EC_POINT z);
void ec_bls509_fp_to_oct(unsigned char *os, size_t *size, const EC_POINT z);
void ec_bls509_fp_from_oct(EC_POINT z, const unsigned char *os, size_t size);
void cat_int_str(unsigned char *os, size_t *oslen, const mpz_t i, const unsigned char *s, const size_t slen);
void ec_bls509_fp_decompose_scalar_init(mpz_t a1, mpz_t a2, mpz_t b1, mpz_t b2, const mpz_t n, const mpz_t l);
void ec_bls509_fp_init_ec_data(EC_GROUP ec);
void ec_bls509_fp_clear_ec_data(EC_GROUP ec);

// //----------------------------------------------
// // declaration function of elliptic curve
// //----------------------------------------------
// void ec_bls509_fp4_point_init(EC_POINT p);
// void ec_bls509_fp4_point_clear(EC_POINT p);
// void ec_bls509_fp4_point_set(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_point_set_str(EC_POINT z, const char* s);
// void ec_bls509_fp4_point_set_infinity(EC_POINT z);
// void ec_bls509_fp4_point_get_str(char *s, const EC_POINT z);
// void ec_bls509_fp4_add(EC_POINT z, const EC_POINT x, const EC_POINT y);
// void ec_bls509_fp4_dob(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_add_formul(EC_POINT z, const EC_POINT x, const EC_POINT y);
// void ec_bls509_fp4_dob_formul(EC_POINT z, const EC_POINT x);
// //void ec_bls509_fp4_add_formul_homo(EC_POINT z, const EC_POINT x, const EC_POINT y);
// //void ec_bls509_fp4_dob_formul_homo(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_neg(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_sub(EC_POINT z, const EC_POINT x, const EC_POINT y);
// void ec_bls509_fp4_mul(EC_POINT z, const mpz_t s, const EC_POINT x);
// //void ec_bls509_fp4_mul_homo(EC_POINT z, const mpz_t s, const EC_POINT x);
// void ec_bls509_fp4_mul_naf(EC_POINT z, const mpz_t s, const EC_POINT x);
// void ec_bls509_fp4_mul_end(EC_POINT z, const mpz_t s, const EC_POINT x);
// void ec_bls509_fp4_frob_p(EC_POINT Q, const EC_POINT P);
// int  ec_bls509_fp4_is_infinity(const EC_POINT P);
// int  ec_bls509_fp4_is_on_curve(const EC_POINT P);
// int  ec_bls509_fp4_cmp(const EC_POINT x, const EC_POINT y);
// void ec_bls509_fp4_make_affine(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_make_affine_homogeneous(EC_POINT z, const EC_POINT x);
// void ec_bls509_fp4_map_to_point(EC_POINT z, const char *s, size_t slen, int t);
// void ec_bls509_fp4_random(EC_POINT z);
// void ec_bls509_fp4_to_oct(unsigned char *os, size_t *size, const EC_POINT z);
// void ec_bls509_fp4_from_oct(EC_POINT z, const unsigned char *os, size_t size);
// void ec_bls509_tw_frob(EC_POINT Q, const EC_POINT P);
// void ec_bls509_tw_frob2(EC_POINT Q, const EC_POINT P);
// void ec_bls509_tw_frob3(EC_POINT Q, const EC_POINT P);
// void ec_bls509_fp4_init_ec_data_aranha(EC_GROUP ec);
// void ec_bls509_fp4_init_ec_data_beuchat(EC_GROUP ec);
// void ec_bls509_fp4_clear_ec_data(EC_GROUP ec);
//
// //----------------------------------------------
// // declaration function of pairing
// //----------------------------------------------
// void ec_bls509_pairing_precomp_beuchat(EC_PAIRING p);
// void ec_bls509_pairing_dob_beuchat(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P);
// void ec_bls509_pairing_add_beuchat(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P, const EC_POINT Q);
// void ec_bls509_pairing_miller_beuchat(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_pairing_finalexp(Element z, const Element x, const EC_PAIRING p);
// void ec_bls509_pairing_beuchat(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_double_pairing_beuchat(Element z, const EC_POINT Q1, const EC_POINT P1, const EC_POINT Q2, const EC_POINT P2, const EC_PAIRING p);
//
// void ec_bls509_pairing_precomp_aranha(EC_PAIRING p);
// void ec_bls509_pairing_dob_aranha_jac(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P);
// void ec_bls509_pairing_add_aranha_jac(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P, const EC_POINT Q);
// void ec_bls509_pairing_dob_aranha_proj(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P);
// void ec_bls509_pairing_add_aranha_proj(EC_POINT T, Element l0, Element l3, Element l4, const EC_POINT P, const EC_POINT Q);
// void ec_bls509_pairing_miller_aranha_jac(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_pairing_miller_aranha_proj(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_pairing_finalexp(Element z, const Element x, const EC_PAIRING p);
// void ec_bls509_pairing_aranha_jac(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_pairing_aranha_proj(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p);
// void ec_bls509_double_pairing_aranha_jac(Element z, const EC_POINT Q1, const EC_POINT P1, const EC_POINT Q2, const EC_POINT P2, const EC_PAIRING p);
// void ec_bls509_double_pairing_aranha_proj(Element z, const EC_POINT Q1, const EC_POINT P1, const EC_POINT Q2, const EC_POINT P2, const EC_PAIRING p);
