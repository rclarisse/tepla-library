//==============================================================
//  extension field ( bls509_fp2 ) implementation with GMP
//--------------------------------------------------------------
//  bls509_fp2 Fp2 := Fp[x]/(x^2 + 1) CLN
//--------------------------------------------------------------
//  2015.10.31 created by kanbara
//==============================================================

#include "ec_bls509_lcl.h"

#define rep(x) (*((mpz_t *)x->data))
#define rep0(x) (((Element *)x->data)[0])
#define rep1(x) (((Element *)x->data)[1])

#define field(x) (x->field)
#define order(x) (x->field->order)
#define beta(x) ((x->field->irre_poly)[0])

#define MODOPT

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bls509_fp2_init(Element x)
{
    x->data = (void *)malloc(sizeof(Element) * 2);

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in fp2 init\n");
        exit(100);
    }

    element_init(rep0(x), field(x)->base);
    element_init(rep1(x), field(x)->base);
}

void bls509_fp2_clear(Element x)
{
    if (x->data != NULL)
    {
        element_clear(rep0(x));
        element_clear(rep1(x));

        free(x->data);
        x->data = NULL;
    }
}

void bls509_fp2_set(Element x, const Element y)
{
    bls509_fp_set(rep0(x), rep0(y));
    bls509_fp_set(rep1(x), rep1(y));
}

void bls509_fp2_set_fp(Element z, const Element x, const Element y)
{
    bls509_fp_set(rep0(z), x);
    bls509_fp_set(rep1(z), y);
}

void bls509_fp2_set_str(Element x, const char *s)
{
    int i = 0;
    int len = strlen(s);

    char msg[260], *p, *c = NULL;

    if (len > 260) {
        fprintf(stderr, "error: input string is too long, string must be smaller than 140\n");
        exit(200);
    }

    strcpy(msg, s);

    p = msg;

    while ((*p) != '\0') {
        if ((*p) == ' ') {
            if (i == 0) {
                c = p;
            }
            i++;
        }
        p++;
    }

    if (i != 1) {
        fprintf(stderr, "error: input string is not correct\n");
        exit(200);
    }

    (*c) = '\0';

    bls509_fp_set_str(rep0(x), msg);
    bls509_fp_set_str(rep1(x), ++c);
}

void bls509_fp2_get_str(char *s, const Element x)
{
    char s1[130], s2[130];

    bls509_fp_get_str(s1, rep0(x));
    bls509_fp_get_str(s2, rep1(x));

    sprintf(s, "%s %s", s1, s2);
}

void bls509_fp2_set_zero(Element x)
{
    bls509_fp_set_zero(rep0(x));
    bls509_fp_set_zero(rep1(x));
}

void bls509_fp2_set_one(Element x)
{
    bls509_fp_set_one(rep0(x));
    bls509_fp_set_zero(rep1(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bls509_fp2_add(Element z, const Element x, const Element y)
{
    bls509_fp_add(rep0(z), rep0(x), rep0(y));
    bls509_fp_add(rep1(z), rep1(x), rep1(y));
}

void bls509_fp2_add_one(Element z, const Element x)
{
    bls509_fp_add_one(rep0(z), rep0(x));
    bls509_fp_set(rep1(z), rep1(x));
}

void bls509_fp2_addn(Element z, const Element x, const Element y)
{
    bls509_fp_addn(rep0(z), rep0(x), rep0(y));
    bls509_fp_addn(rep1(z), rep1(x), rep1(y));
}

void bls509_fp2_neg(Element z, const Element x)
{
    bls509_fp_neg(rep0(z), rep0(x));
    bls509_fp_neg(rep1(z), rep1(x));
}

void bls509_fp2_sub(Element z, const Element x, const Element y)
{
    bls509_fp_sub(rep0(z), rep0(x), rep0(y));
    bls509_fp_sub(rep1(z), rep1(x), rep1(y));
}

void bls509_fp2_subn(Element z, const Element x, const Element y)
{
    bls509_fp_subn(rep0(z), rep0(x), rep0(y));
    bls509_fp_subn(rep1(z), rep1(x), rep1(y));
}

//------------------------------------------------------
//  multipulication is implemented by Karatsuba method.
//------------------------------------------------------
void bls509_fp2_mul(Element z, const Element x, const Element y)
{
    Element* t = field(z)->base->tmp;

    bls509_fp_addn(t[1], rep0(x), rep1(x)); // t1 = x0 + x1
    bls509_fp_addn(t[2], rep0(y), rep1(y)); // t2 = y0 + y1
    bls509_fp_muln(t[0], t[1], t[2]);       // t0 = t1 * t2
    bls509_fp_muln(t[1], rep0(x), rep0(y)); // t1 = x0 * y0
    bls509_fp_muln(t[2], rep1(x), rep1(y)); // t2 = x1 * y1
    bls509_fp_subn(t[0], t[0], t[1]);       //
    bls509_fp_subn(rep1(z), t[0], t[2]);    // z1 = x0*y1 + y0*x1
    bls509_fp_mod(rep1(z), rep1(z));        //
    bls509_fp_subn(rep0(z), t[1], t[2]);    //
    bls509_fp_mod(rep0(z), rep0(z));        // z0 = t1 - t2

}

void bls509_fp2_muln(Element z, const Element x, const Element y)
{
    Element* t = field(z)->base->tmp;

    bls509_fp_muln(t[0], rep0(x), rep0(y)); // t0 = x0 * y0
    bls509_fp_muln(t[1], rep1(x), rep1(y)); // t1 = x0 * y1
    bls509_fp_addn(t[2], rep0(x), rep1(x)); // t2 = x0 + x1
    bls509_fp_addn(t[3], rep0(y), rep1(y)); // t3 = y0 + y1
    bls509_fp_muln(t[4], t[2], t[3]);	      // t4 = t2 * t3
    bls509_fp_addn(t[5], t[0], t[1]);	      // t5 = t0 + t1
    bls509_fp_subn(rep1(z), t[4], t[5]);	  // z1 = t4 - t5
    bls509_fp_sub(t[6], t[0], t[1]);		    // t6 = t0 - t1
    bls509_fp_OP2(rep0(z), t[6]);           // z0 = t6

}


//--------------------------------------------------------
//  multiplication of element of fp and element of fp^2
//--------------------------------------------------------
void bls509_fp2_mul_p(Element z, const Element x, const Element y)
{
    if (field(x)->ID == bls509_fp)
    {
        bls509_fp_mul(rep0(z), x, rep0(y));
        bls509_fp_mul(rep1(z), x, rep1(y));
    }
    else if (field(y)->ID == bls509_fp)
    {
        bls509_fp_mul(rep0(z), rep0(x), y);
        bls509_fp_mul(rep1(z), rep1(x), y);
    }
    else
    {
        fprintf(stderr, "error: input should be element in bls509_fp2\n");
        exit(200);
    }
}

void bls509_fp2_mul_c(Element z, const Element x, const mpz_t c)
{
    bls509_fp_mulc(rep0(z), rep0(x), c);
    bls509_fp_mulc(rep1(z), rep1(x), c);
}

void bls509_fp2_div_2(Element z, const Element x)
{
    bls509_fp_div2(rep0(z), rep0(x));
    bls509_fp_div2(rep1(z), rep1(x));
}

void bls509_fp2_inv(Element z, const Element x)
{
    Element* t = field(z)->base->tmp;

    bls509_fp_muln(t[0], rep0(x), rep0(x)); // t0 = x0^2
    bls509_fp_muln(t[1], rep1(x), rep1(x)); // t1 = x1^2
    bls509_fp_addn(t[0], t[0], t[1]);       // t0 = x0^2 + x1^2
    bls509_fp_inv(t[1], t[0]);              // t1 = (x0^2 + x1^2)^-1
    bls509_fp_mul(rep0(z), rep0(x), t[1]);  // z0 = x0*(x0^2 + x1^2)^-1
    bls509_fp_mul(rep1(z), rep1(x), t[1]);  // z1 = x1*(x0^2 + x1^2)^-1
    bls509_fp_neg(rep1(z), rep1(z));        // z1 = -x1*(x0^2 + x1^2)^-1

}

void bls509_fp2_dob(Element z, const Element x)
{
    bls509_fp_dob(rep0(z), rep0(x));
    bls509_fp_dob(rep1(z), rep1(x));
}

void bls509_fp2_tri(Element z, const Element x)
{
    bls509_fp_tri(rep0(z), rep0(x));
    bls509_fp_tri(rep1(z), rep1(x));
}

// multiplication by the root / indeterminate
void bls509_fp2_xi_mul(Element z, const Element x)
{
    bls509_fp_neg(rep0(z), rep1(x));
    bls509_fp_set(rep1(z), rep0(x));
}

// WHAT THE F*CK IS IT DOING ???
// void bls509_fp2_xi_mul_inv(Element z, const Element x)
// {
//     Element* t = field(z)->base->tmp;
//
//     bls509_fp_add(t[0], rep0(x), rep1(x));
//     bls509_fp_sub(t[1], rep1(x), rep0(x));
//
//     bls509_fp_set(rep0(z), t[0]);
//     bls509_fp_set(rep1(z), t[1]);
// }

void bls509_fp2_sqr(Element z, const Element x)
{
    Element* t = field(z)->base->tmp;

    bls509_fp_addn(t[0], rep1(x), rep1(x)); // t0 = 2*x1
    bls509_fp_muln(t[0], t[0], rep0(x));    // t0 = 2*x1*x0
    bls509_fp_addn(t[1], rep0(x), rep1(x)); // t1 = x0+x1
    bls509_fp_subn(t[2], rep0(x), rep1(x)); // t2 = x0-x1
    bls509_fp_muln(t[1], t[1], t[2]);       // t1 = t1*t2
    bls509_fp_mod(rep1(z), t[0]);           // z1 = t0
    bls509_fp_mod(rep0(z), t[1]);		        // z0 = t1

}

void bls509_fp2_sqrn(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bls509_fp_addn(t[0], rep0(x), rep1(x));
    bls509_fp_sub(t[1], rep0(x), rep1(x));
    bls509_fp_muln(rep0(z), t[0], t[1]);
    bls509_fp_addn(t[0], rep0(x), rep0(x));
    bls509_fp_muln(rep1(z), t[0], rep1(x));
}

void bls509_fp2_pow(Element z, const Element x, const mpz_t exp)
{
    long t, i;
    Element c;

    element_init(c, field(z));
    element_set(c, x);

    t = (long)mpz_sizeinbase(exp, 2);

    for (i = t - 2; i >= 0; i--)
    {
        element_sqr(c, c);
        if (mpz_tstbit(exp, i)) {
            element_mul(c, c, x);
        }
    }

    element_set(z, c);
    element_clear(c);
}

void bls509_fp2_mod(Element z, const Element x)
{
    bls509_fp_mod(rep0(z), rep0(x));
    bls509_fp_mod(rep1(z), rep1(x));
}

void bls509_fp2_OP1_1(Element z, const Element x)
{
    bls509_fp_OP1_1(rep0(z), rep0(x));
    bls509_fp_OP1_1(rep1(z), rep1(x));
}
void bls509_fp2_OP1_2(Element z, const Element x)
{
    bls509_fp_OP1_2(rep0(z), rep0(x));
    bls509_fp_OP1_2(rep1(z), rep1(x));
}
void bls509_fp2_OP2(Element z, const Element x)
{
    bls509_fp_OP2(rep0(z), rep0(x));
    bls509_fp_OP2(rep1(z), rep1(x));
}

void bls509_fp2_frob_p(Element z, const Element x)
{
    bls509_fp_set(rep0(z), rep0(x));
    bls509_fp_neg(rep1(z), rep1(x));
}

void bls509_fp2_conj(Element z, const Element x)
{
    bls509_fp_set(rep0(z), rep0(x));
    bls509_fp_neg(rep1(z), rep1(x));
}

//---------------------------------------------------------
//  precomputation for Fp2 operation
//---------------------------------------------------------
// void bls509_fp2_precomp(Field f)
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

//--------------------------------------------------
//  square root in extended Fp
//  Generalized Atkin algorithm for card = 9 mod 16
//  (see Siguna Müller: on the computation of square roots in finite fields)
//--------------------------------------------------
int bls509_fp2_sqrt(Element z, const Element x)
{
    mpz_t exp, deux;

    Element *t = field(x)->tmp;

    if (!element_is_sqr(x)) {
        return FALSE;
    }

    mpz_init_set_ui(deux, 2);

    bls509_fp2_mul_c(t[0], x, deux); // t0 = 2*x

    mpz_init_set(exp, order(x));

    mpz_sub_ui(exp, exp, 1);
    mpz_fdiv_q_2exp(exp, exp, 2); // exp = (p^2 - 1) / 4

    bls509_fp2_pow(t[1], t[0], exp); // t1 = t0^exp (t1 = +/- 1 mod p^2)

    bls509_fp2_random(t[2]);
    if (bls509_fp2_is_one(t[1])) {
      while (bls509_fp2_is_sqr(t[2]) == TRUE) {
        bls509_fp2_random(t[2]);
      }
    } else {
      while (bls509_fp2_is_sqr(t[2]) == FALSE) {
        bls509_fp2_random(t[2]);
      }
    }

    bls509_fp2_sqr(t[1], t[2]); // t1 = t2^2
    bls509_fp2_mul(t[3], t[1], t[0]); // t3 = t1*t0 = 2*x*t2^2

    mpz_set(exp, order(x));
    mpz_sub_ui(exp, exp, 9);
    mpz_fdiv_q_2exp(exp, exp, 4); // exp = (p^2 - 9) / 16

    bls509_fp2_pow(t[4], t[3], exp); // t4 = t3^exp
    bls509_fp2_mul(t[5], t[3], t[4]);
    bls509_fp2_mul(t[5], t[5], t[4]); // t5 = t3*t4^2 (t5^2 = -1)

    bls509_fp2_mul(z, x, t[4]);
    bls509_fp2_mul(z, z, t[2]);
    bls509_fp2_set_str(t[6], "1 0"); // t6 = 1
    bls509_fp2_sub(t[5], t[5], t[6]);
    bls509_fp2_mul(z, z, t[5]); // z = t4 * t2 * x * (t5 - 1)

    mpz_clear(exp);
    mpz_clear(deux);

    return TRUE;
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bls509_fp2_is_zero(const Element x)
{
    return (bls509_fp_is_zero(rep1(x)) && bls509_fp_is_zero(rep0(x)));
}

int bls509_fp2_is_one(const Element x)
{
    if (bls509_fp_is_zero(rep1(x)))
    {
        return (bls509_fp_is_one(rep0(x)));
    }
    return FALSE;
}

int bls509_fp2_cmp(const Element x, const Element y)
{
    if (bls509_fp_cmp(rep1(x), rep1(y)) == 0)
    {
        if (bls509_fp_cmp(rep0(x), rep0(y)) == 0) {
            return 0;
        }
    }
    return 1;
}

int bls509_fp2_is_sqr(const Element x)
{
    int hr = FALSE;

    Element *t = field(x)->base->tmp;

    if (element_is_zero(x)) {
        return FALSE;
    }

    bls509_fp_inv(t[0], rep1(x));
    bls509_fp_mul(t[0], t[0], rep0(x));
    bls509_fp_sqr(t[0], t[0]);
    bls509_fp_add(t[0], t[0], field(x)->irre_poly[0]);

    hr = bls509_fp_is_sqr(t[0]);

    return hr;
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bls509_fp2_random(Element z)
{
    element_random(rep0(z));
    element_random(rep1(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bls509_fp2_to_mpz(mpz_t a, const Element x)
{
    mpz_set(a, rep(rep0(x)));   // a = rep0
    mpz_addmul(a, rep(rep1(x)), field(x)->base->order);   //a = a + rep1*p
}

void bls509_fp2_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t s0;

    unsigned char b0[128];
    mpz_t z;

    mpz_init(z);

    bls509_fp2_to_mpz(z, x);
    mpz_export(b0, &s0, 1, sizeof(*b0), 1, 0, z);

    memset(os, 0x00, 128);

    memcpy(&os[128 - (int)s0], b0, s0);

    (*size) = 128;

    mpz_clear(z);
}

void bls509_fp2_from_oct(Element x, const unsigned char *os, const size_t size)
{
    mpz_t quo, rem;

    if (size < 128) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_init(quo);
    mpz_init(rem);

    mpz_import(quo, size, 1, sizeof(*os), 1, 0, os);

    mpz_tdiv_qr(quo, rem, quo, field(x)->base->order);
    mpz_set(rep(rep0(x)), rem);
    mpz_set(rep(rep1(x)), quo);

    mpz_clear(quo);
    mpz_clear(rem);
}
