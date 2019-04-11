//==============================================================
//  extension field ( bls509_fp4 ) implementation with GMP
//--------------------------------------------------------------
//  bls509_fp4 Fp4 := Fp2[y]/(y^2 + x + 1)     CNL
//--------------------------------------------------------------
//  2019.04.11 created by rclarisse
//==============================================================

#include "ec_bls509_lcl.h"

#define rep(x) (*((mpz_t *)x->data))
#define rep0(x) (((Element *)x->data)[0])
#define rep1(x) (((Element *)x->data)[1])

#define field(x) (x->field)
#define order(x) (x->field->order)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bls509_fp4_init(Element x)
{
    x->data = (void *)malloc(sizeof(Element) * 2);

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in bls509_fp4 init\n");
        exit(100);
    }

    element_init(rep0(x), field(x)->base);
    element_init(rep1(x), field(x)->base);
}

void bls509_fp4_clear(Element x)
{
    if (x->data != NULL)
    {
        element_clear(rep0(x));
        element_clear(rep1(x));

        free(x->data);
        x->data = NULL;
    }
}

void bls509_fp4_set(Element x, const Element y)
{
    bls509_fp2_set(rep0(x), rep0(y));
    bls509_fp2_set(rep1(x), rep1(y));
}

void bls509_fp4_set_fp2(Element z, const Element x, const Element y)
{
    bls509_fp2_set(rep0(z), x);
    bls509_fp2_set(rep1(z), y);
}

void bls509_fp4_set_str(Element x, const char *s)
{
    int i = 0;
    int len = strlen(s);

    char msg[520], *p, *c[3];

    if (len > 520) {
        fprintf(stderr, "error: input string is too long, string must be smaller than 520\n");
        exit(200);
    }

    strcpy(msg, s);

    p = msg;

    while ((*p) != '\0')
    {
        if ((*p) == ' ') {
            if (i < 3) {
                c[i] = p;
            }
            i++;
        }
        p++;
    }

    if (i != 3) {
        fprintf(stderr, "error: input string is not correct\n");
        exit(200);
    }

    (*c[0]) = '\0';
    (*c[1]) = '\0';
    (*c[2]) = '\0';

    bls509_fp_set_str(rep0(rep0(x)), msg);
    bls509_fp_set_str(rep0(rep1(x)), ++c[0]);
    bls509_fp_set_str(rep1(rep0(x)), ++c[2]);
    bls509_fp_set_str(rep1(rep1(x)), ++c[3]);
}

void bls509_fp4_get_str(char *s, const Element x)
{
    char s0[130], s1[130], s2[130], s3[130];

    bls509_fp_get_str(s0, rep0(rep0(x)));
    bls509_fp_get_str(s1, rep0(rep1(x)));
    bls509_fp_get_str(s2, rep1(rep0(x)));
    bls509_fp_get_str(s3, rep1(rep1(x)));

    sprintf(s, "%s %s %s %s", s0, s1, s2, s3);
}

void bls509_fp4_set_zero(Element x)
{
    bls509_fp2_set_zero(rep0(x));
    bls509_fp2_set_zero(rep1(x));
}

void bls509_fp4_set_one(Element x)
{
    bls509_fp2_set_one(rep0(x));
    bls509_fp2_set_zero(rep1(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bls509_fp4_add(Element z, const Element x, const Element y)
{
    bls509_fp2_add(rep0(z), rep0(x), rep0(y));
    bls509_fp2_add(rep1(z), rep1(x), rep1(y));
}

void bls509_fp4_addn(Element z, const Element x, const Element y)
{
    bls509_fp2_addn(rep0(z), rep0(x), rep0(y));
    bls509_fp2_addn(rep1(z), rep1(x), rep1(y));
}

void bls509_fp4_neg(Element z, const Element x)
{
    bls509_fp2_neg(rep0(z), rep0(x));
    bls509_fp2_neg(rep1(z), rep1(x));
}

void bls509_fp4_sub(Element z, const Element x, const Element y)
{
    bls509_fp2_sub(rep0(z), rep0(x), rep0(y));
    bls509_fp2_sub(rep1(z), rep1(x), rep1(y));
}

void bls509_fp4_subn(Element z, const Element x, const Element y)
{
    bls509_fp2_subn(rep0(z), rep0(x), rep0(y));
    bls509_fp2_subn(rep1(z), rep1(x), rep1(y));
}

// xi is the element x + 1 in Fp2 \sub Fp4
void bls509_fp4_mul(Element z, const Element x, const Element y)
{
    Element *t = field(z)->base->tmp;

    bls509_fp2_mul(t[0], rep0(x), rep0(y));  // t0 = x0*y0
    bls509_fp2_mul(t[1], rep1(x), rep1(y));  // t1 = x1*y1

    //------------------------------------------
    //  z0 = x0*y0 - x1*y1*xi
    //------------------------------------------
    bls509_fp2_xi_mul(t[2], t[1]);       // t2 = x1*y1*xi
    bls509_fp2_sub(rep0(z), t[0], t[2]); // z0 = t0 - t2

    //------------------------------------------
    //  z1 = x0*y1 + x1*y0
    //------------------------------------------
    bls509_fp2_add(t[2], t[0], t[1]);       // t2 = t0 + t1 = x0*y0 + x1*y1
    bls509_fp2_add(t[3], rep0(x), rep1(x)); // t3 = x0 + x1
    bls509_fp2_add(t[4], rep0(y), rep1(y)); // t4 = y0 + y1
    bls509_fp2_mul(t[5], t[3], t[4]);       // t5 = t3*t4 = (x0 + x1)(y0 + y1)
    bls509_fp2_sub(rep1(z), t[5], t[2]);    // z1 = t5 - t2

}

// void bls509_fp4_muln(Element z, const Element x, const Element y)
// {
//     Element *t = field(z)->base->tmp;
//
//     bls509_fp2_muln(t[0], rep0(x), rep0(y));  // t0 = a0*b0
//     bls509_fp2_muln(t[1], rep1(x), rep1(y));  // t1 = a1*b1
//     bls509_fp2_muln(t[2], rep2(x), rep2(y));  // t2 = a2*b2
//     bls509_fp2_OP1_2(t[0], t[0]);
//     bls509_fp2_OP1_2(t[1], t[1]);
//     bls509_fp2_OP1_2(t[2], t[2]);
//     bls509_fp2_addn(t[8], rep1(x), rep2(x));	 // t8 = a1+a2
//     bls509_fp2_addn(t[9], rep1(y), rep2(y));	 // t9 = b1+b2
//     bls509_fp2_muln(t[3], t[8], t[9]);		 // t3 = t8*t9
//     bls509_fp2_OP2(t[3], t[3]);
//     bls509_fp2_addn(t[4], t[1], t[2]);		 // t4 = t1+t2
//     bls509_fp_sub(rep0(t[3]), rep0(t[3]), rep0(t[4]));
//     bls509_fp_OP2(rep0(t[3]), rep0(t[3]));
//     bls509_fp_subn(rep1(t[3]), rep1(t[3]), rep1(t[4]));
//     bls509_fp2_xi_mul(t[4], t[3]);			 // t4 = xi*t3
//     bls509_fp2_OP2(t[4], t[4]);
//     bls509_fp2_add(rep0(z), t[4], t[0]);		 // t5 = t4+t0
//     bls509_fp2_OP2(rep0(z), rep0(z));
//
//     bls509_fp2_addn(t[8], rep0(x), rep1(x));	 // t8 = a0+a1
//     bls509_fp2_addn(t[9], rep0(y), rep1(y));	 // t9 = b0+b1
//     bls509_fp2_muln(t[3], t[8], t[9]);		 // t3 = t8*t9
//     bls509_fp2_OP2(t[3], t[3]);
//     bls509_fp2_addn(t[4], t[0], t[1]);		 // t4 = t0+t1
//     bls509_fp_sub(rep0(t[3]), rep0(t[3]), rep0(t[4]));
//     bls509_fp_OP2(rep0(t[3]), rep0(t[3]));
//     bls509_fp_subn(rep1(t[3]), rep1(t[3]), rep1(t[4]));
//     bls509_fp_sub(rep0(t[4]), rep0(t[2]), rep1(t[2]));
//     bls509_fp_OP1_1(rep0(t[4]), rep0(t[4]));
//     bls509_fp_addn(rep1(t[4]), rep0(t[2]), rep1(t[2]));
//     bls509_fp2_add(rep1(z), t[3], t[4]);		 // t6 = t3+t4
//     bls509_fp2_OP2(rep1(z), rep1(z));
//
//     bls509_fp2_addn(t[8], rep0(x), rep2(x));	 // t8 = a0+a2
//     bls509_fp2_addn(t[9], rep0(y), rep2(y));	 // t9 = b0+b2
//     bls509_fp2_muln(t[3], t[8], t[9]);		 // t3 = t8*t9
//     bls509_fp2_OP2(t[3], t[3]);
//     bls509_fp2_addn(t[4], t[0], t[2]);		 // t4 = t0+t1
//     bls509_fp_sub(rep0(t[3]), rep0(t[3]), rep0(t[4]));
//     bls509_fp_OP2(rep0(t[3]), rep0(t[3]));
//     bls509_fp_subn(rep1(t[3]), rep1(t[3]), rep1(t[4]));
//     bls509_fp_add(rep0(rep2(z)), rep0(t[3]), rep0(t[1]));
//     bls509_fp_OP2(rep0(rep2(z)), rep0(rep2(z)));
//     bls509_fp_addn(rep1(rep2(z)), rep1(t[3]), rep1(t[1]));
// }

//--------------------------------------------------------
//   z = x * gamma ( Fp12 : fp4[x]/x^2-gamma )
//--------------------------------------------------------
// void bls509_fp4_gm_mul(Element z, const Element x)
// {
//     if (z == x) {
//         fprintf(stderr, "fail gm mul\n");
//         exit(500);
//     }
//
//     bls509_fp2_xi_mul(rep0(z), rep2(x));
//     bls509_fp2_set(rep1(z), rep0(x));
//     bls509_fp2_set(rep2(z), rep1(x));
// }

//--------------------------------------------------------
//   z = x * (y, 0)
//--------------------------------------------------------
void bls509_fp4_mul_fp2(Element z, const Element x, const Element y)
{
    if (field(y)->ID != bls509_fp2)
    {
        fprintf(stderr, "error: input should be element in bls509_fp4\n");
        exit(200);
    }

    bls509_fp2_mul(rep0(z), rep0(x), y);
    bls509_fp2_mul(rep1(z), rep1(x), y);
}

//--------------------------------------------------------
//   z = x * (y1, y2, 0)
//--------------------------------------------------------
// void bls509_fp4_mul_fp2_2(Element z, const Element x, const Element y1, const Element y2)
// {
//     Element *v = field(z)->base->tmp;
//
//     if (field(y1)->ID != bls509_fp2 || field(y2)->ID != bls509_fp2)
//     {
//         fprintf(stderr, "error: input should be element in bls509_fp4\n");
//         exit(200);
//     }
//
//     bls509_fp2_mul(v[0], rep0(x), y1);
//     bls509_fp2_mul(v[1], rep1(x), y2);
//     bls509_fp2_add(v[2], rep1(x), rep2(x));
//     bls509_fp2_mul(v[2], v[2], y2);
//     bls509_fp2_sub(v[2], v[2], v[1]);
//     bls509_fp2_xi_mul(v[2], v[2]);
//     bls509_fp2_add(v[2], v[2], v[0]);
//     bls509_fp2_add(rep1(z), rep0(x), rep1(x));
//     bls509_fp2_add(v[3], y1, y2);
//     bls509_fp2_mul(rep1(z), rep1(z), v[3]);
//     bls509_fp2_sub(rep1(z), rep1(z), v[0]);
//     bls509_fp2_sub(rep1(z), rep1(z), v[1]);
//     bls509_fp2_add(rep2(z), rep0(x), rep2(x));
//     bls509_fp2_mul(rep2(z), rep2(z), y1);
//     bls509_fp2_sub(rep2(z), rep2(z), v[0]);
//     bls509_fp2_add(rep2(z), rep2(z), v[1]);
//     bls509_fp2_set(rep0(z), v[2]);
// }

//--------------------------------------------------------
//   z = x * (0, y, 0)
//--------------------------------------------------------
// void bls509_fp4_mul_fp2_3(Element z, const Element x, const Element y)
// {
//     if (field(y)->ID != bls509_fp2)
//     {
//         fprintf(stderr, "error: input should be element in bls509_fp4\n");
//         exit(200);
//     }
//
//     bls509_fp2_mul(rep1(z), rep0(x), y);
//     bls509_fp2_mul(rep2(z), rep1(x), y);
//     bls509_fp2_mul(rep0(z), rep2(x), y);
//     bls509_fp2_xi_mul(rep0(z), rep0(z));
// }

//--------------------------------------------------------
//   z = x * (y1, 0, y2)
//--------------------------------------------------------
// void bls509_fp4_mul_fp2_4(Element z, const Element x, const Element y1, const Element y2)
// {
//     Element *v = field(z)->base->tmp;
//
//     if (field(y1)->ID != bls509_fp2 || field(y2)->ID != bls509_fp2)
//     {
//         fprintf(stderr, "error: input should be element in bls509_fp4\n");
//         exit(200);
//     }
//
//     bls509_fp2_mul(v[0], rep0(x), y1);
//     bls509_fp2_mul(v[1], rep2(x), y2);
//
//     bls509_fp2_add(v[3], rep0(x), rep2(x));
//     bls509_fp2_add(v[4], y1, y2);
//     bls509_fp2_mul(rep2(z), v[3], v[4]);
//     bls509_fp2_sub(rep2(z), rep2(z), v[0]);
//     bls509_fp2_sub(rep2(z), rep2(z), v[1]);
//
//     bls509_fp2_mul(v[3], rep1(x), y2);
//     bls509_fp2_mul(v[4], rep1(x), y1);
//
//     bls509_fp2_xi_mul(v[1], v[1]);
//     bls509_fp2_xi_mul(v[3], v[3]);
//
//     bls509_fp2_add(rep0(z), v[0], v[3]);
//     bls509_fp2_add(rep1(z), v[1], v[4]);
// }

void bls509_fp4_inv(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bls509_fp2_sqr(t[0], rep0(x));    // t0 = x0^2
    bls509_fp2_sqr(t[1], rep1(x));    // t1 = x1^2
    bls509_fp2_xi_mul(t[1], t[1]);    // t1 = x1^2 * xi
    bls509_fp2_add(t[0], t[0], t[1]); // t0 = x0^2 + x1^2 * xi
    bls509_fp2_inv(t[0], t[0]);

    bls509_fp2_mul(rep0(z), rep0(x), t[0]);
    bls509_fp2_mul(rep1(z), rep1(z), t[0]);
    bls509_fp2_neg(rep1(z), rep1(z));
}

void bls509_fp4_mod(Element z, const Element x)
{
    bls509_fp2_mod(rep0(z), rep0(x));
    bls509_fp2_mod(rep1(z), rep1(x));
    bls509_fp2_mod(rep2(z), rep2(x));
}

// void bls509_fp4_OP1_1(Element z, const Element x)
// {
//     bls509_fp2_OP1_1(rep0(z), rep0(x));
//     bls509_fp2_OP1_1(rep1(z), rep1(x));
//     bls509_fp2_OP1_1(rep2(z), rep2(x));
// }
//
// void bls509_fp4_OP1_2(Element z, const Element x)
// {
//     bls509_fp2_OP1_2(rep0(z), rep0(x));
//     bls509_fp2_OP1_2(rep1(z), rep1(x));
//     bls509_fp2_OP1_2(rep2(z), rep2(x));
// }
//
// void bls509_fp4_OP2(Element z, const Element x)
// {
//     bls509_fp2_OP2(rep0(z), rep0(x));
//     bls509_fp2_OP2(rep1(z), rep1(x));
//     bls509_fp2_OP2(rep2(z), rep2(x));
// }

void bls509_fp4_dob(Element z, const Element x)
{
    bls509_fp2_dob(rep0(z), rep0(x));
    bls509_fp2_dob(rep1(z), rep1(x));
}

void bls509_fp4_tri(Element z, const Element x)
{
    bls509_fp2_tri(rep0(z), rep0(x));
    bls509_fp2_tri(rep1(z), rep1(x));
}

void bls509_fp4_sqr(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bls509_fp2_addn(t[0], rep1(x), rep1(x)); // t0 = 2*x1
    bls509_fp2_muln(t[0], t[0], rep0(x));    // t0 = 2*x1*x0
    bls509_fp2_mod(rep1(z), t[0]);           // z1 = t0

    bls509_fp2_sqrn(t[0], rep0(x));    // t0 = x0^2
    bls509_fp2_sqrn(t[1], rep1(x));    // t1 = x1^2
    bls509_fp2_xi_mul(t[1], t[1]);     // t1 = x1^2*xi
    bls509_fp2_subn(t[2], t[0], t[1]); // t2 = x0^2 - x1^2*xi
    bls509_fp2_mod(rep0(z), t[2]);     // z0 = t2
}

/*
void bls509_fp4_frob_p(Element z, const Element x)
{
}
*/

// void bls509_fp4_conj(Element z, const Element x)
// {
//     bls509_fp2_conj(rep0(z), rep0(x));
//     bls509_fp2_conj(rep1(z), rep1(x));
// }

// //---------------------------------------------------------
// //  precomputation for fp4 operation
// //---------------------------------------------------------
// void bls509_fp4_precomp(Field f)
// {
//     field_precomp_p precomp = NULL;
//
//     precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));
//
//     precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
//     bls509_fp2_precomp_sqrt(precomp->ps, f);
//
//     precomp->pf = NULL;
//
//     f->precomp = (void *)precomp;
// }
//
// //---------------------------------------------------------
// //  precomputation for fp4 operation for pairing_init
// //---------------------------------------------------------
// void bls509_fp4_precomp_for_pairing_init(Field f)
// {
//     field_precomp_p precomp = NULL;
//
//     precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));
//
//     precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
//     bls509_fp2_precomp_sqrt_for_fp4init(precomp->ps, f);
//
//     precomp->pf = NULL;
//
//     f->precomp = (void *)precomp;
// }

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bls509_fp4_is_zero(const Element x)
{
    if (bls509_fp2_is_zero(rep1(x)))
    {
        if (bls509_fp2_is_zero(rep0(x)))
        {
            return TRUE;
        }
    }
    return FALSE;
}

int bls509_fp4_is_one(const Element x)
{
    if (bls509_fp2_is_zero(rep1(x)))
    {
        if (bls509_fp2_is_one(rep0(x)))
        {
            return TRUE;
        }
    }
    return FALSE;
}

int bls509_fp4_cmp(const Element x, const Element y)
{
    if (bls509_fp2_cmp(rep1(x), rep1(y)) == 0)
    {
        if (bls509_fp2_cmp(rep0(x), rep0(y)) == 0) {
            return 0;
        }
    }
    return 1;
}

int bls509_fp4_is_sqr(const Element x)
{
    int k = 1;

    Element *t = field(x)->base->tmp;

    if (element_is_zero(x)) {
        return FALSE;
    }

    k *= bls509_fp2_is_sqr(rep2(x)) ? 1 : -1;

    bls509_fp2_sqr(t[1], rep1(x));
    bls509_fp2_mul(t[2], rep0(x), rep2(x));
    bls509_fp2_sub(t[1], t[1], t[2]);      // t1 = x1^2-x0*x2
    bls509_fp2_mul(t[2], rep0(x), rep1(x));
    bls509_fp2_sqr(t[3], rep2(x));
    bls509_fp2_xi_mul(t[3], t[3]);
    bls509_fp2_sub(t[2], t[2], t[3]);      // t2 = x0*x1-x2^2*xi
    bls509_fp2_inv(t[1], t[1]);
    bls509_fp2_mul(t[1], t[1], t[2]);      // t1 = t2 / t1

    bls509_fp2_inv(t[2], rep2(x));
    bls509_fp2_mul(t[3], t[2], rep1(x));
    bls509_fp2_sub(t[3], t[3], t[1]);
    bls509_fp2_mul(t[3], t[3], t[1]);      // t3 = ((x1/x2)-t1)t1

    bls509_fp2_mul(t[2], t[2], rep0(x));
    bls509_fp2_sub(t[2], t[2], t[3]);      // t2 = (x0/x2)-t3

    k *= bls509_fp2_is_sqr(t[2]) ? 1 : -1;

    return (k == 1);
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bls509_fp4_random(Element z)
{
    bls509_fp2_random(rep0(z));
    bls509_fp2_random(rep1(z));
    bls509_fp2_random(rep2(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bls509_fp4_to_mpz(mpz_t a, const Element x)
{
    mpz_mul(a, rep(rep1(rep2(x))), field(x)->base->base->order);   // a = rep12*p
    mpz_add(a, a, rep(rep1(rep1(x))));   // a = a + rep11
    mpz_mul(a, a, field(x)->base->base->order);   // a = a*p
    mpz_add(a, a, rep(rep1(rep0(x))));   // a = a + rep10
    mpz_mul(a, a, field(x)->base->base->order);   // a = a*p
    mpz_add(a, a, rep(rep0(rep2(x))));   // a = a + rep02
    mpz_mul(a, a, field(x)->base->base->order);   // a = a*p
    mpz_add(a, a, rep(rep0(rep1(x))));   // a = a + rep01
    mpz_mul(a, a, field(x)->base->base->order);   // a = a*p
    mpz_add(a, a, rep(rep0(rep0(x))));   // a = a + rep00
}

void bls509_fp4_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t s0;

    unsigned char b0[190];
    mpz_t z;

    mpz_init(z);

    bls509_fp4_to_mpz(z, x);
    mpz_export(b0, &s0, 1, sizeof(*b0), 1, 0, z);

    memset(os, 0x00, 190);

    memcpy(&os[190 - (int)s0], b0, s0);

    (*size) = 190;

    mpz_clear(z);
}

void bls509_fp4_from_oct(Element x, const unsigned char *os, const size_t size)
{
    mpz_t quo, rem;

    if (size < 190) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_init(quo);
    mpz_init(rem);

    mpz_import(quo, size, 1, sizeof(*os), 1, 0, os);

    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep0(rep0(x))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep0(rep1(x))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep0(rep2(x))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep1(rep0(x))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep1(rep1(x))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->order);
    mpz_set(rep(rep1(rep2(x))), rem);

    mpz_clear(quo);
    mpz_clear(rem);
}
