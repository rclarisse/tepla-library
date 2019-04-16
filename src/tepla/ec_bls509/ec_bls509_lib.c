//==================================================
// Function which making instance of structures
//==================================================
#include "ec_bls509_lcl.h"

#define TMP_NUM 10

#define SAFE_FREE(p) { if(p){ free(p); (p)=NULL; } }

//----------------------------------------------
//  function release field "bls509"
//----------------------------------------------
void ec_bls509_field_clear(Field f)
{
    unsigned int i, j;

    if (f->precomp != NULL)
    {
        field_precomp_p precomp = (field_precomp_p)(f->precomp);

        field_precomp_frob_p pf = precomp->pf;

        if (pf != NULL)
        {
            for (i = 0; i < pf->glen1; i++) {
                element_clear(pf->gamma1[i]);
            }
            for (i = 0; i < pf->glen2; i++) {
                element_clear(pf->gamma2[i]);
            }
            for (i = 0; i < pf->glen3; i++) {
                element_clear(pf->gamma3[i]);
            }
            free(pf->gamma1);
            free(pf->gamma2);
            free(pf->gamma3);
            free(pf);
        }
        SAFE_FREE(f->precomp);
    }

    if (f->irre_poly != NULL)
    {
        for (i = 0, j = f->irre_poly_num; i < j; i++) {
            element_clear(f->irre_poly[i]);
        }
        f->irre_poly_num = 0;
        f->irre_poly_deg = 0;
        SAFE_FREE(f->irre_poly);
    }

    if (f->tmp != NULL)
    {
        for (i = 0; i < TMP_NUM; i++) {
            element_clear(f->tmp[i]);
        }
        SAFE_FREE(f->tmp);
    }

    if (f->base != NULL)
    {
        field_clear(f->base);
        SAFE_FREE(f->base);
    }

    mpz_clear(f->order);
    mpz_clear(f->OP1_1);
    mpz_clear(f->OP1_2);
    mpz_clear(f->OP2);

    SAFE_FREE(f->field_name);

    f->str_len = 0;
    f->oct_len = 0;

    f->init = NULL;
    f->clear = NULL;
    f->set = NULL;
    f->set_str = NULL;
    f->get_str = NULL;
    f->set_zero = NULL;
    f->set_one = NULL;
    f->add = NULL;
    f->sub = NULL;
    f->neg = NULL;
    f->mul = NULL;
    f->sqr = NULL;
    f->inv = NULL;
    f->pow = NULL;
    f->sqrt = NULL;
    f->is_zero = NULL;
    f->is_one = NULL;
    f->is_sqr = NULL;
    f->cmp = NULL;
    f->random = NULL;
    f->to_oct = NULL;
    f->from_oct = NULL;
}

//----------------------------------------------
//  function creating field bls509_fp
//----------------------------------------------
void ec_bls509_fp_new(Field f)
{
    int i;

    f->type = Field_fp;

    set_field_name(f, "bls509_fp");

    f->ID = bls509_fp;

    f->str_len = 130;
    f->oct_len = 64;

    //------------------------------
    //  set pointer of function
    //------------------------------
    f->init     = bls509_fp_init;
    f->clear    = bls509_fp_clear;
    f->set      = bls509_fp_set;
    f->set_str  = bls509_fp_set_str;
    f->get_str  = bls509_fp_get_str;
    f->set_zero = bls509_fp_set_zero;
    f->set_one  = bls509_fp_set_one;

    f->add  = bls509_fp_add;
    f->sub  = bls509_fp_sub;
    f->neg  = bls509_fp_neg;
    f->mul  = bls509_fp_mul;
    f->sqr  = bls509_fp_sqr;
    f->inv  = bls509_fp_inv;
    f->pow  = bls509_fp_pow;
    f->sqrt = bls509_fp_sqrt;

    f->is_zero = bls509_fp_is_zero;
    f->is_one  = bls509_fp_is_one;
    f->is_sqr  = bls509_fp_is_sqr;
    f->cmp     = bls509_fp_cmp;

    f->random = bls509_fp_random;

    f->to_oct   = bls509_fp_to_oct;
    f->from_oct = bls509_fp_from_oct;

    //-----------------------------------------
    //  base field
    //-----------------------------------------
    f->base = NULL;

    //-----------------------------------------
    //  characteristic of prime field
    //-----------------------------------------
    //    p = (t - 1)^2 * (t^8 - t^4 + 1) / 3 + t
    //    t = -1 + 2^11 - 2^28 - 2^51 = -2251800082118657
    //-----------------------------------------
    mpz_init_set_str(f->order, "155556FFFF39CA9BFCEDF2B4F9C0ECF6CB8AC8495D187E8C32EA0103E01090BB626E85BF7C18A0F0CFCB5C6071BAD3D2EE63BD076E8D9300A13D118DB8BFD2AB", 16);
    mpz_init_set_str(f->OP1_1, "0", 16);
    mpz_init_set_str(f->OP1_2, "0", 16);
    mpz_init_set_str(f->OP2, "0", 16);

    //-----------------------------------------
    //  irreducible polynomial
    //-----------------------------------------
    f->irre_poly_num = 0;
    f->irre_poly_deg = 1;
    f->irre_poly = NULL;

    //-----------------------------------------
    //  pre-computation for square root
    //-----------------------------------------
    f->precomp = NULL;

    //----------------------------------
    //  temporary element init
    //----------------------------------
    f->tmp = (Element *)malloc(sizeof(Element) * TMP_NUM);
    for (i = 0; i < TMP_NUM; i++) {
        element_init(f->tmp[i], f);
    }

    return;
}

//----------------------------------------------
//  function creating field bls509_fp2
//----------------------------------------------
void ec_bls509_fp2_new(Field f)
{
    int i;

    f->type = Field_fpn;

    set_field_name(f, "bls509_fp2");

    f->ID = bls509_fp2;

    f->str_len = 260;
    f->oct_len = 128;

    //------------------------------
    //  set pointer of function
    //------------------------------
    f->init     = bls509_fp2_init;
    f->clear    = bls509_fp2_clear;
    f->set      = bls509_fp2_set;
    f->set_str  = bls509_fp2_set_str;
    f->get_str  = bls509_fp2_get_str;
    f->set_zero = bls509_fp2_set_zero;
    f->set_one  = bls509_fp2_set_one;

    f->add  = bls509_fp2_add;
    f->sub  = bls509_fp2_sub;
    f->neg  = bls509_fp2_neg;
    f->mul  = bls509_fp2_mul;
    f->sqr  = bls509_fp2_sqr;
    f->inv  = bls509_fp2_inv;
    f->pow  = bls509_fp2_pow;
    f->sqrt = bls509_fp2_sqrt;

    f->is_zero = bls509_fp2_is_zero;
    f->is_one  = bls509_fp2_is_one;
    f->is_sqr  = bls509_fp2_is_sqr;
    f->cmp     = bls509_fp2_cmp;

    f->random = bls509_fp2_random;

    f->to_oct   = bls509_fp2_to_oct;
    f->from_oct = bls509_fp2_from_oct;

    //-----------------------------------------
    //  set base field
    //-----------------------------------------
    f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));

    field_init(f->base, "bls509_fp");

    //-----------------------------------------
    //  characteristic of prime field
    //-----------------------------------------
    //    p = (t - 1)^2 * (t^8 - t^4 + 1) / 3 + t
    //    t = -1 + 2^11 - 2^28 - 2^51 = -2251800082118657
    //-----------------------------------------
    mpz_init(f->order);
    mpz_mul(f->order, f->base->order, f->base->order);

    mpz_init_set_str(f->OP1_1, "0", 16);
    mpz_init_set_str(f->OP1_2, "0", 16);
    mpz_init_set_str(f->OP2, "0", 16);

    //-----------------------------------------
    //  pre-computation for square root
    //-----------------------------------------
    f->precomp = NULL;

    //-----------------------------------------
    //  Irreducible polynomial: x^2 + 1
    //-----------------------------------------
    f->irre_poly_num = 1;
    f->irre_poly_deg = 2;

    f->irre_poly = (Element *)malloc(sizeof(Element));
    element_init(f->irre_poly[0], f->base);
    element_set_str(f->irre_poly[0], "1");

    //----------------------------------
    //  temporary element init
    //----------------------------------
    f->tmp = (Element *)malloc(sizeof(Element) * TMP_NUM);
    for (i = 0; i < TMP_NUM; i++) {
        element_init(f->tmp[i], f);
    }

    return;
}

//----------------------------------------------
//  function creating field bls509_fp4
//----------------------------------------------
void ec_bls509_fp4_new(Field f)
{
    int i;

    f->type = Field_fpn;

    set_field_name(f, "bls509_fp4");

    f->ID = bls509_fp4;

    f->str_len = 520;
    f->oct_len = 256;

    //------------------------------
    //  set pointer of function
    //------------------------------
    f->init     = bls509_fp4_init;
    f->clear    = bls509_fp4_clear;
    f->set      = bls509_fp4_set;
    f->set_str  = bls509_fp4_set_str;
    f->get_str  = bls509_fp4_get_str;
    f->set_zero = bls509_fp4_set_zero;
    f->set_one  = bls509_fp4_set_one;

    f->add  = bls509_fp4_add;
    f->sub  = bls509_fp4_sub;
    f->neg  = bls509_fp4_neg;
    f->mul  = bls509_fp4_mul;
    f->sqr  = bls509_fp4_sqr;
    f->inv  = bls509_fp4_inv;
    f->pow  = bls509_fp4_pow;
    f->sqrt = NULL; // bls509_fp4_sqrt;

    f->is_zero = bls509_fp4_is_zero;
    f->is_one  = bls509_fp4_is_one;
    f->is_sqr  = bls509_fp4_is_sqr;

    f->cmp = bls509_fp4_cmp;

    f->random = bls509_fp4_random;

    f->to_oct   = bls509_fp4_to_oct;
    f->from_oct = bls509_fp4_from_oct;

    //-----------------------------------------
    //  set base field
    //-----------------------------------------
    f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));

    field_init(f->base, "bls509_fp2");

    //-----------------------------------------
    //  characteristic of prime field
    //-----------------------------------------
    //    p = (t - 1)^2 * (t^8 - t^4 + 1) / 3 + t
    //    t = -1 + 2^11 - 2^28 - 2^51 = -2251800082118657
    //-----------------------------------------
    mpz_init(f->order);
    mpz_mul(f->order, f->base->order, f->base->order);

    mpz_init_set_str(f->OP1_1, "0", 16);
    mpz_init_set_str(f->OP1_2, "0", 16);
    mpz_init_set_str(f->OP2, "0", 16);

    //-----------------------------------------
    //  pre-computation for square root
    //-----------------------------------------
    f->precomp = NULL;

    //-----------------------------------------
    //  irreducible polynomial: y^2 + x + 1
    //-----------------------------------------
    f->irre_poly_num = 1;
    f->irre_poly_deg = 2;

    f->irre_poly = (Element *)malloc(sizeof(Element));

    element_init(f->irre_poly[0], f->base);
    element_set_str(f->irre_poly[0], "1 1");

    //----------------------------------
    //  temporary element init
    //----------------------------------
    f->tmp = (Element *)malloc(sizeof(Element) * TMP_NUM);
    for (i = 0; i < TMP_NUM; i++) {
        element_init(f->tmp[i], f);
    }

    return;
}

//----------------------------------------------
//  function creating field bls509_fp12
//----------------------------------------------
// void ec_bls509_fp12_new(Field f)
// {
//     int i;
//
//     f->type = Field_fpn;
//
//     set_field_name(f, "bls509_fp12");
//
//     f->ID = bls509_fp12;
//
//     f->str_len = 1560;
//     f->oct_len = 768;
//
//     //------------------------------
//     //  set pointer of function
//     //------------------------------
//     f->init     = bls509_fp12_init;
//     f->clear    = bls509_fp12_clear;
//     f->set      = bls509_fp12_set;
//     f->set_str  = bls509_fp12_set_str;
//     f->get_str  = bls509_fp12_get_str;
//     f->set_zero = bls509_fp12_set_zero;
//     f->set_one  = bls509_fp12_set_one;
//
//     f->add  = bls509_fp12_add;
//     f->sub  = bls509_fp12_sub;
//     f->neg  = bls509_fp12_neg;
//     f->mul  = bls509_fp12_mul;
//     f->sqr  = bls509_fp12_sqr;
//     f->inv  = bls509_fp12_inv;
//     f->pow  = bls509_fp12_pow_naf;
//     f->sqrt = bls509_fp12_sqrt;
//
//     f->is_zero = bls509_fp12_is_zero;
//     f->is_one  = bls509_fp12_is_one;
//     f->is_sqr  = bls509_fp12_is_sqr;
//     f->cmp     = bls509_fp12_cmp;
//
//     f->random = bls509_fp12_random;
//
//     f->to_oct   = bls509_fp12_to_oct;
//     f->from_oct = bls509_fp12_from_oct;
//
//     //-----------------------------------------
//     //  set base field
//     //-----------------------------------------
//     f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
//     field_init(f->base, "bls509_fp4");
//
//     //-----------------------------------------
//     //  characteristic of prime field
//     //-----------------------------------------
//     //    p = (t - 1)^2 * (t^8 - t^4 + 1) / 3 + t
//     //    t = -1 + 2^11 - 2^28 - 2^51 = -2251800082118657
//     //-----------------------------------------
//     mpz_init(f->order);
//     mpz_mul(f->order, f->base->order, f->base->order);
//     mpz_mul(f->order, f->order, f->base->order);
//
//     mpz_init_set_str(f->OP1_1, "0", 16);
//     mpz_init_set_str(f->OP1_2, "0", 16);
//     mpz_init_set_str(f->OP2, "0", 16);
//
//     //-----------------------------------------
//     //  irreducible polynomial: z^3 + y
//     //-----------------------------------------
//     f->irre_poly_num = 1;
//     f->irre_poly_deg = 3;
//
//     f->irre_poly = (Element *)malloc(sizeof(Element));
//
//     element_init(f->irre_poly[0], f->base);
//     element_set_str(f->irre_poly[0], "0 1 0 0");
//
//     //----------------------------------
//     //  temporary element init
//     //----------------------------------
//     f->tmp = (Element *)malloc(sizeof(Element) * TMP_NUM);
//     for (i = 0; i < TMP_NUM; i++) {
//         element_init(f->tmp[i], f);
//     }
//
//     return;
// }


//----------------------------------------------
//  function creating field bls509_fp24
//----------------------------------------------
// void ec_bls509_fp24_new(Field f)
// {
//     int i;
//
//     f->type = Field_fpn;
//
//     set_field_name(f, "bls509_fp24");
//
//     f->ID = bls509_fp24;
//
//     f->str_len = 3120;
//     f->oct_len = 1536;
//
//     //------------------------------
//     //  set pointer of function
//     //------------------------------
//     f->init     = bls509_fp24_init;
//     f->clear    = bls509_fp24_clear;
//     f->set      = bls509_fp24_set;
//     f->set_str  = bls509_fp24_set_str;
//     f->get_str  = bls509_fp24_get_str;
//     f->set_zero = bls509_fp24_set_zero;
//     f->set_one  = bls509_fp24_set_one;
//
//     f->add  = bls509_fp24_add;
//     f->sub  = bls509_fp24_sub;
//     f->neg  = bls509_fp24_neg;
//     f->mul  = bls509_fp24_mul;
//     f->sqr  = bls509_fp24_sqr;
//     f->inv  = bls509_fp24_inv;
//     f->pow  = bls509_fp24_pow;
//     f->sqrt = bls509_fp24_sqrt;
//
//     f->is_zero = bls509_fp24_is_zero;
//     f->is_one  = bls509_fp24_is_one;
//     f->is_sqr  = bls509_fp24_is_sqr;
//
//     f->cmp = bls509_fp24_cmp;
//
//     f->random = bls509_fp24_random;
//
//     f->to_oct   = bls509_fp24_to_oct;
//     f->from_oct = bls509_fp24_from_oct;
//
//     //-----------------------------------------
//     //  set base field
//     //-----------------------------------------
//     f->base = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
//     field_init(f->base, "bls509_fp12");
//
//     //-----------------------------------------
//     //  characteristic of prime field
//     //-----------------------------------------
//     //    p = (t - 1)^2 * (t^8 - t^4 + 1) / 3 + t
//     //    t = -1 + 2^11 - 2^28 - 2^51 = -2251800082118657
//     //-----------------------------------------
//     mpz_init(f->order);
//     mpz_mul(f->order, f->base->order, f->base->order);
//
//     mpz_init_set_str(f->OP1_1, "0", 16);
//     mpz_init_set_str(f->OP1_2, "0", 16);
//     mpz_init_set_str(f->OP2, "0", 16);
//
//     //-----------------------------------------
//     //  irreducible polynomial: y^2 + x
//     //-----------------------------------------
//     f->irre_poly_num = 1;
//     f->irre_poly_deg = 2;
//
//     f->irre_poly = (Element *)malloc(sizeof(Element));
//
//     element_init(f->irre_poly[0], f->base);
//     element_set_str(f->irre_poly[0], "0 1");
//
//     //----------------------------------
//     //  temporary element init
//     //----------------------------------
//     f->tmp = (Element *)malloc(sizeof(Element) * TMP_NUM);
//     for (i = 0; i < TMP_NUM; i++) {
//         element_init(f->tmp[i], f);
//     }
//
//     return;
// }

//----------------------------------------------
//  function generating elliptic curve method
//----------------------------------------------
// void ec_bls509_fp_method_new(EC_METHOD method)
// {
//     method->point_init  = ec_bls509_fp_point_init;
//     method->point_clear = ec_bls509_fp_point_clear;
//     method->point_set = ec_bls509_fp_point_set;
//     method->point_set_infinity = ec_bls509_fp_point_set_infinity;
//     method->point_set_str = ec_bls509_fp_point_set_str;
//     method->point_get_str = ec_bls509_fp_point_get_str;
//
//     method->add = ec_bls509_fp_add;
//     method->dob = ec_bls509_fp_dob;
//     method->neg = ec_bls509_fp_neg;
//     method->sub = ec_bls509_fp_sub;
//
// #ifdef ENABLE_FASTALG
//     method->mul = ec_bls509_fp_mul_end;
// #else
//     method->mul = ec_bls509_fp_mul_naf;
// #endif
//
//     method->is_infinity = ec_bls509_fp_is_infinity;
//     method->is_on_curve = ec_bls509_fp_is_on_curve;
//     method->cmp = ec_bls509_fp_cmp;
//
//     method->make_affine = ec_bls509_fp_make_affine;
//     method->map_to_point = ec_bls509_fp_map_to_point;
//     method->random = ec_bls509_fp_random;
//     method->to_oct = ec_bls509_fp_to_oct;
//     method->from_oct = ec_bls509_fp_from_oct;
// }


//----------------------------------------------
//  function generating elliptic curve group
//----------------------------------------------
// void ec_bls509_fp_group_new(EC_GROUP ec)
// {
//     ec->type = Curve_BLS;
//
//     set_curve_name(ec, "ec_bls509_fp");
//
//     ec->ID = ec_bls509_fp;
//
//     ec->str_len = 132;
//     ec->oct_len = 65;
//
//     ec->field = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
//     field_init(ec->field, "bls509_fp");
//
//     ec->method = (struct ec_method_st *)malloc(sizeof(struct ec_method_st));
//     ec_bls509_fp_method_new(ec->method);
//
//     point_init(ec->generator, ec);
//     point_set_str(ec->generator, "[1,d45589b158faaf6ab0e4ad38d998e9982e7ff63964ee1460342a592677cccb0]");
//
//     element_init(ec->a, ec->field);
//     element_init(ec->b, ec->field);
//
//     element_set_zero(ec->a);
//     element_set_str(ec->b, "5");
//
//     mpz_init_set_str(ec->order, "2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001", 16);
//     mpz_init_set_str(ec->trace, "5F408FD0060000000000000000000001", 16);
//     mpz_init_set_str(ec->cofactor, "1", 16);
//
//     ec_bls509_fp_init_ec_data(ec);
// }

//----------------------------------------------
//  function generating elliptic curve method
//----------------------------------------------
// void ec_bls509_tw_method_new(EC_METHOD method)
// {
//     method->point_init  = ec_bls509_fp_point_init;
//     method->point_clear = ec_bls509_fp_point_clear;
//     method->point_set = ec_bls509_fp_point_set;
//     method->point_set_infinity = ec_bls509_fp_point_set_infinity;
//     method->point_set_str = ec_bls509_fp4_point_set_str;
//     method->point_get_str = ec_bls509_fp4_point_get_str;
//
//     method->add = ec_bls509_fp4_add;
//     method->dob = ec_bls509_fp4_dob;
//     method->neg = ec_bls509_fp4_neg;
//     method->sub = ec_bls509_fp4_sub;
//
// #ifdef ENABLE_FASTALG
//     method->mul = ec_bls509_fp4_mul_end;
// #else
//     method->mul = ec_bls509_fp4_mul;
// #endif
//
//     method->is_infinity = ec_bls509_fp_is_infinity;
//     method->is_on_curve = ec_bls509_fp4_is_on_curve;
//     method->cmp = ec_bls509_fp_cmp;
//
//     method->make_affine = ec_bls509_fp4_make_affine;
//     method->map_to_point = ec_bls509_fp4_map_to_point;
//     method->random = ec_bls509_fp4_random;
//     method->to_oct = ec_bls509_fp4_to_oct;
//     method->from_oct = ec_bls509_fp4_from_oct;
// }

//----------------------------------------------
//  function generating elliptic curve group
//----------------------------------------------
// void ec_bls509_tw_group_new(EC_GROUP ec)
// {
//     ec->type = Curve_BN;
//
//     set_curve_name(ec, "ec_bls509_tw");
//
//     ec->ID = ec_bls509_fp4;
//
//     ec->str_len = 262;
//     ec->oct_len = 129;
//
//     ec->field = (struct ec_field_st *)malloc(sizeof(struct ec_field_st));
//     field_init(ec->field, "bls509_fp4");
//
//     ec->method = (struct ec_method_st *)malloc(sizeof(struct ec_method_st));
//     ec_bls509_tw_method_new(ec->method);
//
//     point_init(ec->generator, ec);
//     point_set_str(ec->generator, "[19b0bea4afe4c330da93cc3533da38a9f430b471c6f8a536e81962ed967909b5 a1cf585585a61c6e9880b1f2a5c539f7d906fff238fa6341e1de1a2e45c3f72,17abd366ebbd65333e49c711a80a0cf6d24adf1b9b3990eedcc91731384d2627 0ee97d6de9902a27d00e952232a78700863bc9aa9be960C32f5bf9fd0a32d345]");
//
//     element_init(ec->a, ec->field);
//     element_init(ec->b, ec->field);
//
//     element_set_zero(ec->a);
//     element_set_str(ec->b, "0 -1");
//
//     mpz_init_set_str(ec->order, "2370FB049D410FBE4E761A9886E502411DC1AF70120000017E80600000000001", 16);
//     mpz_init_set_str(ec->trace, "2370FB049D410FBDC02400000000000000000000000000000000000000000001", 16);
//     mpz_init_set_str(ec->cofactor, "2370FB049D410FBE4E761A9886E50241DC42CF101E0000017E80600000000001", 16);
//
//     //ec_bls509_fp4_init_ec_data(ec);
// }

//----------------------------------------------
//  clear curve group : ec_bls509
//----------------------------------------------
// void ec_bls509_group_clear(EC_GROUP ec)
// {
//     point_clear(ec->generator);
//
//     element_clear(ec->a);
//     element_clear(ec->b);
//
//     mpz_clear(ec->order);
//     mpz_clear(ec->trace);
//     mpz_clear(ec->cofactor);
//
//     if (ec->ID == ec_bls509_fp)
//     {
//         ec_bls509_fp_clear_ec_data(ec);
//     }
//     else if (ec->ID == ec_bls509_fp4)
//     {
//         ec_bls509_fp4_clear_ec_data(ec);
//     }
//
//     SAFE_FREE(ec->method);
//
//     field_clear(ec->field);
//     SAFE_FREE(ec->field);
//
//     SAFE_FREE(ec->curve_name);
//
//     ec->str_len = 0;
//     ec->oct_len = 0;
// }

//-------------------------------------------
// pairing group : Init, Clear
//-------------------------------------------
// void ec_bls509_pairing_new(EC_PAIRING p)
// {
//     p->type = Pairing_ECBLS509;
//
//     set_pairing_name(p, "ECBLS509");
//
//     p->pairing = ec_bls509_pairing_beuchat;
//     p->pairing_double = ec_bls509_double_pairing_beuchat;
//
//     curve_init(p->g1, "ec_bls509_fp");
//     curve_init(p->g2, "ec_bls509_tw");
//
//     //ec_bls509_fp24_new_for_pairing_init(p->g3);
//     p->g3->field_init  = ec_bls509_fp24_new;
//     p->g3->field_clear = ec_bls509_field_clear;
//
//     //ec_bls509_pairing_precomp(p);
// }

// void ec_bls509_pairing_clear(EC_PAIRING p)
// {
//     if (p->precomp != NULL)
//     {
//         pairing_precomp_p precomp = (pairing_precomp_p)(p->precomp);
//
//         free(precomp->si);
//         free(precomp->ti);
//         free(precomp);
//         p->precomp = NULL;
//     }
//
//     p->pairing = NULL;
//
//     curve_clear(p->g1);
//     curve_clear(p->g2);
//     field_clear(p->g3);
// }
