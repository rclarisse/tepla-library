//==============================================================
//  finite field ( bls509_fp24 ) implementation with GMP
//--------------------------------------------------------------
//  bls509_fp24 Fp24 := Fp12[z]/(t^2 + z) CNL
//--------------------------------------------------------------
//  2019.04.19 created by rclarisse
//==============================================================

#include "ec_bls509_lcl.h"

#define rep(x) (*((mpz_t *)x->data))
#define rep0(x) (((Element *)x->data)[0])
#define rep1(x) (((Element *)x->data)[1])
#define rep2(x) (((Element *)x->data)[2])

#define field(x) (x->field)
#define order(x) (x->field->order)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bls509_fp24_init(Element x)
{
    x->data = (void *)malloc(sizeof(Element) * 2);

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in bls509_fp24 init\n");
        exit(100);
    }

    element_init(rep0(x), field(x)->base);
    element_init(rep1(x), field(x)->base);
}

void bls509_fp24_clear(Element x)
{
    if (x->data != NULL)
    {
        element_clear(rep0(x));
        element_clear(rep1(x));

        free(x->data);
        x->data = NULL;
    }
}

void bls509_fp24_set(Element x, const Element y)
{
    bls509_fp12_set(rep0(x), rep0(y));
    bls509_fp12_set(rep1(x), rep1(y));
}

void bls509_fp24_set_fp12(Element z, const Element x, const Element y)
{
    bls509_fp12_set(rep0(z), x);
    bls509_fp12_set(rep1(z), y);
}

void bls509_fp24_set_str(Element x, const char *s)
{
    int i = 0;
    int len = strlen(s);

    char msg[3120], *p, *c[23];

    if (len > 3120) {
        fprintf(stderr, "error: input string is too long, string must be smaller than 3120\n");
        exit(200);
    }

    strcpy(msg, s);

    p = msg;

    while ((*p) != '\0')
    {
        if ((*p) == ' ') {
            if (i < 23) {
                c[i] = p;
            }
            i++;
        }
        p++;
    }

    if (i != 23) {
        fprintf(stderr, "error: input string is not correct %d\n", i);
        exit(200);
    }

    (*c[0]) = '\0';
    (*c[1]) = '\0';
    (*c[2]) = '\0';
    (*c[3]) = '\0';
    (*c[4]) = '\0';
    (*c[5]) = '\0';
    (*c[6]) = '\0';
    (*c[7]) = '\0';
    (*c[8]) = '\0';
    (*c[9]) = '\0';
    (*c[10]) = '\0';
    (*c[11]) = '\0';
    (*c[12]) = '\0';
    (*c[13]) = '\0';
    (*c[14]) = '\0';
    (*c[15]) = '\0';
    (*c[16]) = '\0';
    (*c[17]) = '\0';
    (*c[18]) = '\0';
    (*c[19]) = '\0';
    (*c[20]) = '\0';
    (*c[21]) = '\0';
    (*c[22]) = '\0';

    bls509_fp_set_str(rep0(rep0(rep0(rep0(x)))), msg);
    bls509_fp_set_str(rep1(rep0(rep0(rep0(x)))), ++c[0]);
    bls509_fp_set_str(rep0(rep1(rep0(rep0(x)))), ++c[1]);
    bls509_fp_set_str(rep1(rep1(rep0(rep0(x)))), ++c[2]);
    bls509_fp_set_str(rep0(rep0(rep1(rep0(x)))), ++c[3]);
    bls509_fp_set_str(rep1(rep0(rep1(rep0(x)))), ++c[4]);
    bls509_fp_set_str(rep0(rep1(rep1(rep0(x)))), ++c[5]);
    bls509_fp_set_str(rep1(rep1(rep1(rep0(x)))), ++c[6]);
    bls509_fp_set_str(rep0(rep0(rep2(rep0(x)))), ++c[7]);
    bls509_fp_set_str(rep1(rep0(rep2(rep0(x)))), ++c[8]);
    bls509_fp_set_str(rep0(rep1(rep2(rep0(x)))), ++c[9]);
    bls509_fp_set_str(rep1(rep1(rep2(rep0(x)))), ++c[10]);
    bls509_fp_set_str(rep0(rep0(rep0(rep1(x)))), ++c[11]);
    bls509_fp_set_str(rep1(rep0(rep0(rep1(x)))), ++c[12]);
    bls509_fp_set_str(rep0(rep1(rep0(rep1(x)))), ++c[13]);
    bls509_fp_set_str(rep1(rep1(rep0(rep1(x)))), ++c[14]);
    bls509_fp_set_str(rep0(rep0(rep1(rep1(x)))), ++c[15]);
    bls509_fp_set_str(rep1(rep0(rep1(rep1(x)))), ++c[16]);
    bls509_fp_set_str(rep0(rep1(rep1(rep1(x)))), ++c[17]);
    bls509_fp_set_str(rep1(rep1(rep1(rep1(x)))), ++c[18]);
    bls509_fp_set_str(rep0(rep0(rep2(rep1(x)))), ++c[19]);
    bls509_fp_set_str(rep1(rep0(rep2(rep1(x)))), ++c[20]);
    bls509_fp_set_str(rep0(rep1(rep2(rep1(x)))), ++c[21]);
    bls509_fp_set_str(rep1(rep1(rep2(rep1(x)))), ++c[22]);
}

void bls509_fp24_get_str(char *s, const Element x)
{
    char s0[130], s1[130], s2[130], s3[130], s4[130], s5[130], s6[130], s7[130], s8[130], s9[130], s10[130], s11[130];
    char s12[130], s13[130], s14[130], s15[130], s16[130], s17[130], s18[130], s19[130], s20[130], s21[130];
    char s22[130], s23[130];

    bls509_fp_get_str(s0,  rep0(rep0(rep0(rep0(x)))));
    bls509_fp_get_str(s1,  rep1(rep0(rep0(rep0(x)))));
    bls509_fp_get_str(s2,  rep0(rep1(rep0(rep0(x)))));
    bls509_fp_get_str(s3,  rep1(rep1(rep0(rep0(x)))));
    bls509_fp_get_str(s4,  rep0(rep0(rep1(rep0(x)))));
    bls509_fp_get_str(s5,  rep1(rep0(rep1(rep0(x)))));
    bls509_fp_get_str(s6,  rep0(rep1(rep1(rep0(x)))));
    bls509_fp_get_str(s7,  rep1(rep1(rep1(rep0(x)))));
    bls509_fp_get_str(s8,  rep0(rep0(rep2(rep0(x)))));
    bls509_fp_get_str(s9,  rep1(rep0(rep2(rep0(x)))));
    bls509_fp_get_str(s10, rep0(rep1(rep2(rep0(x)))));
    bls509_fp_get_str(s11, rep1(rep1(rep2(rep0(x)))));
    bls509_fp_get_str(s12, rep0(rep0(rep0(rep1(x)))));
    bls509_fp_get_str(s13, rep1(rep0(rep0(rep1(x)))));
    bls509_fp_get_str(s14, rep0(rep1(rep0(rep1(x)))));
    bls509_fp_get_str(s15, rep1(rep1(rep0(rep1(x)))));
    bls509_fp_get_str(s16, rep0(rep0(rep1(rep1(x)))));
    bls509_fp_get_str(s17, rep1(rep0(rep1(rep1(x)))));
    bls509_fp_get_str(s18, rep0(rep1(rep1(rep1(x)))));
    bls509_fp_get_str(s19, rep1(rep1(rep1(rep1(x)))));
    bls509_fp_get_str(s20, rep0(rep0(rep2(rep1(x)))));
    bls509_fp_get_str(s21, rep1(rep0(rep2(rep1(x)))));
    bls509_fp_get_str(s22, rep0(rep1(rep2(rep1(x)))));
    bls509_fp_get_str(s23, rep1(rep1(rep2(rep1(x)))));

    sprintf(s, "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18, s19, s20, s21, s22, s23);
}

void bls509_fp24_set_zero(Element x)
{
    bls509_fp12_set_zero(rep0(x));
    bls509_fp12_set_zero(rep1(x));
}

void bls509_fp24_set_one(Element x)
{
    bls509_fp12_set_one(rep0(x));
    bls509_fp12_set_zero(rep1(x));
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bls509_fp24_add(Element z, const Element x, const Element y)
{
    bls509_fp12_add(rep0(z), rep0(x), rep0(y));
    bls509_fp12_add(rep1(z), rep1(x), rep1(y));
}

void bls509_fp24_neg(Element z, const Element x)
{
    bls509_fp12_neg(rep0(z), rep0(x));
    bls509_fp12_neg(rep1(z), rep1(x));
}

void bls509_fp24_sub(Element z, const Element x, const Element y)
{
    bls509_fp12_sub(rep0(z), rep0(x), rep0(y));
    bls509_fp12_sub(rep1(z), rep1(x), rep1(y));
}

void bls509_fp24_mul(Element z, const Element x, const Element y)
{
    Element *t = field(z)->base->tmp;

    bls509_fp12_mul(t[0], rep0(x), rep0(y));  // t0 = x0*y0
    bls509_fp12_mul(t[1], rep1(x), rep1(y));  // t1 = x1*y1

    bls509_fp12_add(t[2], rep0(x), rep1(x)); // t2 = x0 + x1
    bls509_fp12_add(t[3], rep0(y), rep1(y)); // t3 = y0 + y1

    //------------------------------------------
    //  z0 = x0*y0 - x1*y1*gamma
    //------------------------------------------
    bls509_fp12_gamma_mul(t[4], t[1]);   // t4 = x1*y1*gamma
    bls509_fp12_sub(rep0(z), t[0], t[4]); // z0 = x0*y0 - x1*y1*gamma

    //------------------------------------------
    //  z1 = x0*y1 + x1*y0
    //------------------------------------------
    bls509_fp12_add(t[4], t[0], t[1]);   // t4 = t0 + t1 = x0*y0 + x1*y1
    bls509_fp12_mul(t[5], t[2], t[3]);   // t5 = t2*t3 = (x0 + x1)*(y0 + y1)
    bls509_fp12_sub(rep1(z), t[5], t[4]); // z1 = t5 - t4
}

//----------------------------------------------------------
//  Input  : z in Fp12 and l0, l3, l4 in Fp2
//  Output : z *= { (x0, 0, 0), (x1, x2, 0) } in Fp12
//----------------------------------------------------------
// void bls509_fp24_mul_L(Element z, Element x0, Element x1, Element x2)
// {
//     Element *v = field(z)->base->tmp;
//
//     bls509_fp6_mul_fp2(v[0], rep0(z), x0);
//     bls509_fp6_mul_fp2_2(v[1], rep1(z), x1, x2);
//     bls509_fp6_add(rep1(z), rep1(z), rep0(z));
//     bls509_fp2_add(x0, x0, x1);
//     bls509_fp6_mul_fp2_2(rep1(z), rep1(z), x0, x2);
//     bls509_fp6_gm_mul(rep0(z), v[1]);
//     bls509_fp6_add(rep0(z), rep0(z), v[0]);
//     bls509_fp6_sub(rep1(z), rep1(z), v[0]);
//     bls509_fp6_sub(rep1(z), rep1(z), v[1]);
// }

//----------------------------------------------------------
//  Input  : z in Fp12 and l0, l2, l4 in Fp2
//  Output : z *= { (x0, 0, x1), (0, x2, 0) } in Fp12
//----------------------------------------------------------
// void bls509_fp24_mul_L2(Element z, Element x0, Element x1, Element x2)
// {
//     Element *v = field(z)->base->tmp;
//
//     bls509_fp6_mul_fp2_4(v[0], rep0(z), x0, x1);
//     bls509_fp6_mul_fp2_3(v[1], rep1(z), x2);
//
//     bls509_fp6_add(rep1(z), rep1(z), rep0(z));
//     bls509_fp6_set_fp2(v[2], x0, x2, x1);
//     bls509_fp6_mul(rep1(z), rep1(z), v[2]);
//
//     bls509_fp6_sub(rep1(z), rep1(z), v[0]);
//     bls509_fp6_sub(rep1(z), rep1(z), v[1]);
//     bls509_fp6_gm_mul(rep0(z), v[1]);
//     bls509_fp6_add(rep0(z), rep0(z), v[0]);
// }

void bls509_fp24_inv(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bls509_fp12_sqr(t[0], rep0(x));    // t0 = x0^2
    bls509_fp12_sqr(t[1], rep1(x));    // t1 = x1^2
    bls509_fp12_gamma_mul(t[2], t[1]); // t2 = x1^2*gamma
    bls509_fp12_add(t[0], t[0], t[2]);  // t0 = x0^2 + x1^2*gamma
    bls509_fp12_inv(t[0], t[0]);       // t0 = (x0^2 + x1^2*gamma)^-1

    bls509_fp12_mul(rep0(z), rep0(x), t[0]);
    bls509_fp12_mul(rep1(z), rep1(x), t[0]);
    bls509_fp12_neg(rep1(z), rep1(z));
}

void bls509_fp24_dob(Element z, const Element x)
{
    bls509_fp12_dob(rep0(z), rep0(x));
    bls509_fp12_dob(rep1(z), rep1(x));
}

void bls509_fp24_tri(Element z, const Element x)
{
    bls509_fp12_tri(rep0(z), rep0(x));
    bls509_fp12_tri(rep1(z), rep1(x));
}

void bls509_fp24_sqr(Element z, const Element x)
{
    Element *t = field(z)->base->tmp;

    bls509_fp12_sqr(t[0], rep0(x)); // t0 = x0^2
    bls509_fp12_sqr(t[1], rep1(x)); // t1 = x1^2

    bls509_fp12_mul(t[2], rep0(x), rep1(x)); // t2 = x0*x1

    //------------------------------------------
    //  z0 = x0^2 - gamma*x1^2
    //------------------------------------------
    bls509_fp12_gamma_mul(t[3], t[1]);
    bls509_fp12_sub(rep0(z), t[0], t[3]);

    //------------------------------------------
    //  z1 = 2*x0*x1
    //------------------------------------------
    bls509_fp12_dob(rep1(z), t[2]);
}

//-----------------------------------------------------------
//  exponentiation z = x^exp
//-----------------------------------------------------------
void bls509_fp24_pow(Element z, const Element x, const mpz_t exp)
{
    long t, i;
    Element c;

    element_init(c, field(z));
    element_set(c, x);

    t = mpz_sizeinbase(exp, 2);

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

//-----------------------------------------------------------
//  exponentiation z = x^exp with NAF
//-----------------------------------------------------------
void bls509_fp24_pow_naf(Element z, const Element x, const mpz_t exp)
{
    long t, i;
    Element c, ix;

    int *naf, nlen;

    element_init(c, field(z));
    element_init(ix, field(z));

    element_set(c, x);
    element_inv(ix, x);

    t = mpz_sizeinbase(exp, 2);

    naf = (int *)malloc(sizeof(int) * (t + 1));

    generate_naf(naf, &nlen, exp);

    for (i = nlen - 2; i >= 0; i--)
    {
        element_sqr(c, c);
        if (naf[i])
        {
            if (naf[i] < 0) {
                element_mul(c, c, ix);
            }
            else {
                element_mul(c, c, x);
            }
        }
    }

    element_set(z, c);
    element_clear(c);
    element_clear(ix);

    free(naf);
}

//-----------------------------------------------------------
//  Frobenius Map in Fp12
//-----------------------------------------------------------
//  frob(x) == x^p
//  z = g + h*w  : g = g0 + g1*v + g2*v^2
//               : h = h0 + h1*v + h2*v^2
//-----------------------------------------------------------
// void bls509_fp24_frob_p(Element z, const Element x)
// {
//     field_precomp_frob_p pf;
//
//     pf = ((field_precomp_p)(field(z)->precomp))->pf;
//
//     bls509_fp6_conj(rep0(z), rep0(x));
//     bls509_fp6_conj(rep1(z), rep1(x));
//
//     if (strcmp(x->field->field_name, "bls509_fp24a") == 0)
//     {
//         bls509_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma1)[0]);   //t2 = t2*gamma1
//         bls509_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma1)[1]);   //t3 = t3*gamma2
//         bls509_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma1)[2]);   //t4 = t4*gamma3
//         bls509_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma1)[3]);   //t5 = t5*gamma4
//         bls509_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma1)[4]);   //t6 = t6*gamma5
//     }
//
//     if (strcmp(x->field->field_name, "bls509_fp24b") == 0)
//     {
//         bls509_fp2_mul(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma1)[0]);   //t2 = t2*gamma1
//         bls509_fp2_mul(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma1)[1]);   //t3 = t3*gamma2
//         bls509_fp2_mul(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma1)[2]);   //t4 = t4*gamma3
//         bls509_fp2_mul(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma1)[3]);   //t5 = t5*gamma4
//         bls509_fp2_mul(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma1)[4]);   //t6 = t6*gamma5
//     }
// }
//
// void bls509_fp24_frob_p2(Element z, const Element x)
// {
//     field_precomp_frob_p pf;
//
//     pf = ((field_precomp_p)(field(z)->precomp))->pf;
//
//     bls509_fp2_set(rep0(rep0(z)), rep0(rep0(x)));   //t1 = g0;
//
//     if (strcmp(x->field->field_name, "bls509_fp24a") == 0)
//     {
//         bls509_fp2_mul_p(rep0(rep1(z)), rep0(rep1(x)), (pf->gamma2)[0]);   //t2 = h0*gamma1
//         bls509_fp2_mul_p(rep1(rep0(z)), rep1(rep0(x)), (pf->gamma2)[1]);   //t3 = g1*gamma2
//         bls509_fp2_mul_p(rep1(rep1(z)), rep1(rep1(x)), (pf->gamma2)[2]);   //t4 = h1*gamma3
//         bls509_fp2_mul_p(rep2(rep0(z)), rep2(rep0(x)), (pf->gamma2)[3]);   //t5 = g2*gamma4
//         bls509_fp2_mul_p(rep2(rep1(z)), rep2(rep1(x)), (pf->gamma2)[4]);   //t6 = h2*gamma5
//     }
//
//     if (strcmp(x->field->field_name, "bls509_fp24b") == 0)
//     {
//         bls509_fp2_mul(rep0(rep1(z)), rep0(rep1(x)), (pf->gamma2)[0]);   //t2 = h0*gamma1
//         bls509_fp2_mul(rep1(rep0(z)), rep1(rep0(x)), (pf->gamma2)[1]);   //t3 = g1*gamma2
//         bls509_fp2_mul(rep1(rep1(z)), rep1(rep1(x)), (pf->gamma2)[2]);   //t4 = h1*gamma3
//         bls509_fp2_mul(rep2(rep0(z)), rep2(rep0(x)), (pf->gamma2)[3]);   //t5 = g2*gamma4
//         bls509_fp2_mul(rep2(rep1(z)), rep2(rep1(x)), (pf->gamma2)[4]);   //t6 = h2*gamma5
//     }
// }
//
// void bls509_fp24_frob_p3(Element z, const Element x)
// {
//     field_precomp_frob_p pf;
//
//     pf = ((field_precomp_p)(field(z)->precomp))->pf;
//
//     bls509_fp6_conj(rep0(z), rep0(x));
//     bls509_fp6_conj(rep1(z), rep1(x));
//
//     if (strcmp(x->field->field_name, "bls509_fp24a") == 0)
//     {
//         bls509_fp2_mul_p(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma3)[0]);   //t2 = t2*gamma1
//         bls509_fp2_mul_p(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma3)[1]);   //t3 = t3*gamma2
//         bls509_fp2_mul_p(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma3)[2]);   //t4 = t4*gamma3
//         bls509_fp2_mul_p(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma3)[3]);   //t5 = t5*gamma4
//         bls509_fp2_mul_p(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma3)[4]);   //t6 = t6*gamma5
//     }
//
//     if (strcmp(x->field->field_name, "bls509_fp24b") == 0)
//     {
//         bls509_fp2_mul(rep0(rep1(z)), rep0(rep1(z)), (pf->gamma3)[0]);   //t2 = t2*gamma1
//         bls509_fp2_mul(rep1(rep0(z)), rep1(rep0(z)), (pf->gamma3)[1]);   //t3 = t3*gamma2
//         bls509_fp2_mul(rep1(rep1(z)), rep1(rep1(z)), (pf->gamma3)[2]);   //t4 = t4*gamma3
//         bls509_fp2_mul(rep2(rep0(z)), rep2(rep0(z)), (pf->gamma3)[3]);   //t5 = t5*gamma4
//         bls509_fp2_mul(rep2(rep1(z)), rep2(rep1(z)), (pf->gamma3)[4]);   //t6 = t6*gamma5
//     }
// }
//
// void bls509_fp24_conj(Element z, const Element x)
// {
//     bls509_fp6_set(rep0(z), rep0(x));
//     bls509_fp6_neg(rep1(z), rep1(x));
// }

//------------------------------------------------------------
//  square operation of Fp^4 (for bls509_sqr_forpairing )
//------------------------------------------------------------
// void bls509_fp12_sqr(Element c0, Element c1, const Element a0, const Element a1)
// {
//     Element *t = field(c0)->tmp;
//
//     bls509_fp2_sqr(t[0], a0);   //t0 = a0^2
//     bls509_fp2_sqr(t[1], a1);   //t1 = a1^2
//
//     bls509_fp2_xi_mul(c0, t[1]);  //c0 = t1*xi
//     bls509_fp2_add(c0, c0, t[0]); //c0 = c0 + t0
//
//     bls509_fp2_add(c1, a0, a1);   //c1 = a0 + a1
//     bls509_fp2_sqr(c1, c1);       //c1 = c1^2
//     bls509_fp2_sub(c1, c1, t[0]); //c1 = c1 - t0
//     bls509_fp2_sub(c1, c1, t[1]); //c1 = c1 - t1
// }

//------------------------------------------------------------------
//  Special exponentiation in Fp12 by Karabina
//    input  : element x
//    output : z = x^t : t = -(2^62+2^55+1)
//-------------------------------------------------------------------
//
// void bls509_fp24_sqr_forpairing_karabina(Element z, const Element x)
// {
//
//     Element *T = field(z)->base->base->tmp;
//     //------------------------------------------------
//     //	g = (g0+g1*s)+(g2+g3*s)*t+(g4+g5*s)t^2
//     //	g^2 = h = (h0+h1*s)+(h2+h3*s)*t+(h4+h5*s)t^2
//     //	C(g) = [g2,g3,g4,g5]
//     //	D(g^2) = (h0+h1*s)+(h2+h3*s)*t+(h4+h5*s)t^2
//     //	if h2 neq 0
//     //	h1 = (h5^2*xi+3*h4^2-2*h3)/4x2, x0 = (2*h1^2+h2h5-h3h4)*xi+1
//     //	if h2 eq 0
//     //	h1 = 2*h4h5/h3, h0 = (2*h1^2-3*h3h4)*xi+1
//     // 	h2 = 2*g2+3(S4,5-S4,S5)*xi
//     // 	h3 = 3(S4+S5*xi)-2*g3
//     //	h4 = 3(S2+S3*xi)-2*g4
//     // 	h5 = 2*g5+3(S2,3-S2-S3)
//     //	Si,j = (gi+gj)^2
//     //	Si = gi^2
//     //------------------------------------------------
//
//     // g2 = rep0(rep1(x))
//     // g3 = rep2(rep0(x))
//     // g4 = rep1(rep0(x))
//     // g5 = rep2(rep1(x))
//
//     bls509_fp2_sqr(T[0], rep1(rep0(x)));						// T0 = g4^2
//     bls509_fp2_sqr(T[1], rep2(rep1(x)));						// T1 = g5^2
//     bls509_fp2_add(T[5], rep1(rep0(x)), rep2(rep1(x)));		// t0 = g4+g5
//     bls509_fp2_sqr(T[2], T[5]);								// T2 = t0^2
//     bls509_fp2_add(T[3], T[0], T[1]);						// T3 = T0+T1
//     bls509_fp2_sub(T[3], T[2], T[3]);						// T3 = T2-T3
//     bls509_fp2_mod(T[5], T[3]);								// t0 = T3 mod p
//     bls509_fp2_add(T[6], rep0(rep1(x)), rep2(rep0(x)));		// t1 = g2+g3
//     bls509_fp2_sqr(T[3], T[6]);								// T3 = t1^2
//     bls509_fp2_sqr(T[2], rep0(rep1(x)));						// T2 = g2^2
//     bls509_fp2_xi_mul(T[6], T[5]);							// t1 = t0*xi
//     bls509_fp2_add(T[5], T[6], rep0(rep1(x)));				// t0 = t1+g2
//     bls509_fp2_dob(T[5], T[5]);								// t0 = 2*t0
//     bls509_fp2_add(rep0(rep1(z)), T[5], T[6]);				// c2 = t0+t1
//
//     bls509_fp2_xi_mul(T[4], T[1]);							// T4 = T1*xi
//     bls509_fp2_add(T[4], T[0], T[4]);						// T4 = T0+T4
//     bls509_fp2_mod(T[5], T[4]);								// t0 = T4 mod p
//     bls509_fp2_sub(T[6], T[5], rep2(rep0(x)));				// t1 = t0-g3
//     bls509_fp2_dob(T[6], T[6]);								// t1 = 2*t1
//
//     bls509_fp2_sqr(T[1], rep2(rep0(x)));						// T1 = g3^2
//     bls509_fp2_add(rep2(rep0(z)), T[6], T[5]);				// c3 = t1+t0
//
//     bls509_fp2_xi_mul(T[4], T[1]);							// T4 = T1*xi
//     bls509_fp2_add(T[4], T[2], T[4]);						// T4 = T2+T4
//     bls509_fp2_mod(T[5], T[4]);								// t0 = T4 mod p
//     bls509_fp2_sub(T[6], T[5], rep1(rep0(x)));				// t1 = t0-g4
//     bls509_fp2_dob(T[6], T[6]);								// t1 = 2*t1
//     bls509_fp2_add(rep1(rep0(z)), T[6], T[5]);				// c4 = t1+t0
//
//     bls509_fp2_addn(T[0], T[2], T[1]);						// T0 = T2+T1
//     bls509_fp2_sub(T[3], T[3], T[0]);						// T3 = T3-T0
//     bls509_fp2_mod(T[5], T[3]);								// t0 = T3 mod p
//     bls509_fp2_add(T[6], T[5], rep2(rep1(x)));				// t1 = t0+g5
//     bls509_fp2_dob(T[6], T[6]);								// t1 = 2*t1
//     bls509_fp2_add(rep2(rep1(z)), T[6], T[5]);				// c5 = t1+t0
// }
//
// void bls509_fp24_decompose_forpairing_karabina(Element z, const Element x)
// {
//     Element *t = field(z)->base->base->tmp;
//
//     // g2 = rep0(rep1(x))
//     // g3 = rep2(rep0(x))
//     // g4 = rep1(rep0(x))
//     // g5 = rep2(rep1(x))
//
//     if (bls509_fp2_is_zero(rep0(rep1(x))))
//     {
//         bls509_fp2_mul(t[1], rep1(rep0(x)), rep2(rep1(x)));	// t1 = g4*g5
//         bls509_fp2_dob(t[1], t[1]);							// t1 = 2g4g5
//         bls509_fp2_inv(t[2], rep2(rep0(x)));					// t2 = 1/g3
//         bls509_fp2_mul(t[7], t[1], t[2]); 					// c1 = 2g4g5/g3
//
//         bls509_fp2_sqr(t[2], t[1]);							// t2 = g1^2
//         bls509_fp2_dob(t[2], t[2]);							// t2 = 2t2
//         bls509_fp2_mul(t[3], rep2(rep0(x)), rep1(rep0(x)));	// t3 = g3*g4
//         bls509_fp2_dob(t[4], t[3]);							// t4 = 2t3
//         bls509_fp2_add(t[4], t[4], t[3]);					// t4 = 3t3
//         bls509_fp2_sub(t[0], t[3], t[4]);					// g0 = t3-t4
//         bls509_fp2_xi_mul(t[0], t[0]);						// t0 = t0*xi
//         bls509_fp2_add_one(t[0], t[0]);						// t0 = t0+1
//     }
//     else
//     {
//         bls509_fp2_sqr(t[5], rep2(rep1(x)));					// t5 = g5^2
//         bls509_fp2_xi_mul(t[5], t[5]);						// t5 = g5^2*xi
//         bls509_fp2_sqr(t[4], rep1(rep0(x)));					// t4 = g4^2
//         bls509_fp2_dob(t[6], t[4]);							// t4 = 2*g4^2
//         bls509_fp2_add(t[4], t[6], t[4]);					// t4 = 3*g4^2
//         bls509_fp2_dob(t[3], rep2(rep0(x)));					// t3 = 2*g3
//         bls509_fp2_dob(t[2], rep0(rep1(x)));					// t2 = 2*g2
//         bls509_fp2_dob(t[2], t[2]);							// t2 = 4*g2
//         bls509_fp2_inv(t[2], t[2]);							// t2 = 1/4*g2
//         bls509_fp2_add(t[1], t[5], t[4]);					// t1 = t5+t4
//         bls509_fp2_sub(t[1], t[1], t[3]);					// t1 = t5+t4-t3
//         bls509_fp2_mul(t[7], t[1], t[2]);					// c1 = (t5+t4+t3)*t2
//
//         bls509_fp2_sqr(t[1], t[7]);							// t1 = t1^2
//         bls509_fp2_dob(t[1], t[1]);							// t1 = 2*t1
//         bls509_fp2_mul(t[2], rep0(rep1(x)), rep2(rep1(x)));	// t2 = g2*g5
//         bls509_fp2_mul(t[3], rep2(rep0(x)), rep1(rep0(x)));	// t3 = g3*g4
//         bls509_fp2_dob(t[4], t[3]);							// t4 = 2*t3
//         bls509_fp2_add(t[3], t[4], t[3]);					// t3 = 3*t3
//         bls509_fp2_add(t[0], t[1], t[2]);					// t0 = t1+t2
//         bls509_fp2_sub(t[0], t[0], t[3]); 					// t0 = t0+t3
//         bls509_fp2_xi_mul(t[0], t[0]);						// t0 = t0*xi
//         bls509_fp2_add_one(t[0], t[0]); 						// c0 = t0+1
//     }
//
//     bls509_fp2_set(rep0(rep0(z)), t[0]);					// z0 = c0
//     bls509_fp2_set(rep1(rep0(z)), rep1(rep0(x)));		// z1 = c4
//     bls509_fp2_set(rep2(rep0(z)), rep2(rep0(x)));		// z2 = c3
//     bls509_fp2_set(rep0(rep1(z)), rep0(rep1(x)));		// z3 = c2
//     bls509_fp2_set(rep1(rep1(z)), t[7]);					// z4 = c1
//     bls509_fp2_set(rep2(rep1(z)), rep2(rep1(x)));		// z5 = c5
//
// }
//
// void bls509_fp24_pow_forpairing_karabina(Element z, const Element x, const int *t, int tlen)
// {
//     int i;
//
//     Element z55, z62;
//     element_init(z55, field(z));
//     element_init(z62, field(z));
//
//     if (z == x)
//     {
//         fprintf(stderr, "error: bad input for bls509_fp24_pow_forpairing\n");
//         exit(200);
//     }
//
//     element_set(z, x);
//     for (i = 0; i < tlen - 1; i++)
//     {
//         bls509_fp24_sqr_forpairing_karabina(z, z);
//         if (i == 55 - 1) {
//             bls509_fp24_decompose_forpairing_karabina(z55, z);
//         }
//     }
//     bls509_fp24_decompose_forpairing_karabina(z62, z);
//
//     bls509_fp24_mul(z, x, z55);
//     bls509_fp24_mul(z, z, z62);
//
//     bls509_fp24_conj(z, z);
//
//     element_clear(z55);
//     element_clear(z62);
// }

//------------------------------------------------------------------
//  Special exponentiation in Fp12 by Beuchat
//    input  : element x
//    output : z = x^t : t = 2^62-2^54+2^44
//-------------------------------------------------------------------
//
// void bls509_fp24_sqr_forpairing_beuchat(Element z, const Element x)
// {
//     Element *t = field(z)->base->base->tmp;
//     Element *c = field(z)->base->tmp;
//
//     //------------------------
//     // z = g + h*w
//     // g = g0 + g1*v + g2*v^2
//     // h = h0 + h1*v + h2*v^2
//     //------------------------
//     bls509_fp12_sqr(t[2], t[3], rep0(rep0(x)), rep1(rep1(x)));   //t00, t11 = (g0 + h1*V)^2
//     bls509_fp12_sqr(t[4], t[5], rep0(rep1(x)), rep2(rep0(x)));   //t01, t12 = (h0 + g2*V)^2
//     bls509_fp12_sqr(t[6], t[7], rep1(rep0(x)), rep2(rep1(x)));   //t02, t10 = (g1 + h2*V)^2
//
//     bls509_fp2_xi_mul(t[7], t[7]);       // t10 = t10*xi
//
//     bls509_fp2_tri(t[2], t[2]);          // c00 = 3*t00
//     bls509_fp2_dob(t[8], rep0(rep0(x))); // tmp = 2*g0
//     bls509_fp2_sub(t[2], t[2], t[8]);    // c00 = c00 - tmp
//
//     bls509_fp2_tri(t[4], t[4]);          // c01 = 3*t01
//     bls509_fp2_dob(t[8], rep1(rep0(x))); // tmp = 2*g1
//     bls509_fp2_sub(t[4], t[4], t[8]);    // c01 = c01 - tmp
//
//     bls509_fp2_tri(t[6], t[6]);          // c02 = 3*t02
//     bls509_fp2_dob(t[8], rep2(rep0(x))); // tmp = 2*g2
//     bls509_fp2_sub(t[6], t[6], t[8]);    // c02 = c02 - tmp
//
//     bls509_fp2_tri(t[7], t[7]);          // c10 = 3*t10
//     bls509_fp2_dob(t[8], rep0(rep1(x))); // tmp = 2*h0
//     bls509_fp2_add(t[7], t[7], t[8]);    // c10 = c10 + tmp
//
//     bls509_fp2_tri(t[3], t[3]);          // c11 = 3*t11
//     bls509_fp2_dob(t[8], rep1(rep1(x))); // tmp = 2*h1
//     bls509_fp2_add(t[3], t[3], t[8]);    // c11 = c11 + tmp
//
//     bls509_fp2_tri(t[5], t[5]);          // c12 = 3*t12
//     bls509_fp2_dob(t[8], rep2(rep1(x))); // tmp = 2*h2
//     bls509_fp2_add(t[5], t[5], t[8]);    // c12 = c12 + tmp
//
//     bls509_fp6_set_fp2(c[0], t[2], t[4], t[6]); // c0 = c00 + c01*v + c02*v^2
//     bls509_fp6_set_fp2(c[1], t[7], t[3], t[5]); // c1 = c10 + c11*v + c12*v^2
//
//     bls509_fp24_set_fp6(z, c[0], c[1]);   // z = c0 + c1*w
// }
//
// void bls509_fp24_pow_forpairing_beuchat(Element z, const Element x, const int *t, int tlen)
// {
//     int i;
//     Element ix;
//
//     if (z == x)
//     {
//         fprintf(stderr, "error: bad input for bls509_fp24_pow_forpairing\n");
//         exit(200);
//     }
//
//     element_init(ix, field(x));
//
//     bls509_fp24_set(z, x);
//     bls509_fp24_conj(ix, x);
//
//     for (i = tlen - 2; i >= 0; i--)
//     {
//         bls509_fp24_sqr_forpairing_beuchat(z, z);
//
//         if (t[i])
//         {
//             if (t[i] < 0) {
//                 bls509_fp24_mul(z, z, ix);
//             }
//             else {
//                 bls509_fp24_mul(z, z, x);
//             }
//         }
//     }
//
//     element_clear(ix);
// }
//
// //---------------------------------------------------------
// //  precomputation for Fp12 frobenius
// //---------------------------------------------------------
// void bls509_fp24_precomp_frob_aranha(field_precomp_frob_p pf, const Field f)
// {
//     int i;
//     mpz_t  t, u;
//     Element xi, tmp1, tmp2;
//
//     //---------------------------------
//     //  allocate buffer
//     //---------------------------------
//     Element *g1 = (Element *)malloc(sizeof(Element) * 5);
//     Element *g2 = (Element *)malloc(sizeof(Element) * 5);
//     Element *g3 = (Element *)malloc(sizeof(Element) * 5);
//
//     struct ec_field_st *fp2 = f->base->base;
//
//     //---------------------------------
//     //  initialization
//     //---------------------------------
//     for (i = 0; i < 5; i++)
//     {
//         element_init(g1[i], fp2);
//         element_init(g2[i], fp2);
//         element_init(g3[i], fp2);
//         element_init(tmp2, fp2);
//     }
//
//     //---------------------------------
//     //  xi = xi^{(p-1)/6}
//     //---------------------------------
//     mpz_init(t);
//     mpz_init(u);
//     mpz_sub_ui(t, fp2->base->order, 1);
//     mpz_fdiv_q_ui(t, t, 6);
//
//     element_init(xi, fp2);
//     element_init(tmp1, fp2);
//     element_set_str(xi, "1 1");
//     element_pow(tmp1, xi, t);
//
//     //---------------------------------
//     //  set gamma 1, 2, 3
//     //---------------------------------
//     element_set(g1[0], tmp1);
//
//     for (i = 1; i < 5; i++) {
//         element_mul(g1[i], g1[i - 1], g1[0]);
//     }
//     for (i = 0; i < 5; i++)
//     {
//         bls509_fp2_conj(tmp2, g1[i]);
//         element_mul(g2[i], g1[i], tmp2);
//     }
//     for (i = 0; i < 5; i++) {
//         element_mul(g3[i], g1[i], g2[i]);
//     }
//
//     pf->gamma1 = g1;
//     pf->gamma2 = g2;
//     pf->gamma3 = g3;
//
//     pf->glen1 = pf->glen2 = pf->glen3 = 5;
//
//     mpz_clear(t);
//     mpz_clear(u);
//
//     element_clear(xi);
//     element_clear(tmp1);
//     element_clear(tmp2);
// }
//
// void bls509_fp24_precomp_frob_beuchat(field_precomp_frob_p pf, const Field f)
// {
//     int i;
//     mpz_t  t, u;
//     Element xi, tmp;
//
//     //---------------------------------
//     //  allocate buffer
//     //---------------------------------
//     Element *g1 = (Element *)malloc(sizeof(Element) * 5);
//     Element *g2 = (Element *)malloc(sizeof(Element) * 5);
//     Element *g3 = (Element *)malloc(sizeof(Element) * 5);
//
//     struct ec_field_st *fp  = f->base->base->base;
//     struct ec_field_st *fp2 = f->base->base;
//
//     //---------------------------------
//     //  initialization
//     //---------------------------------
//     for (i = 0; i < 5; i++)
//     {
//         element_init(g1[i], fp);
//         element_init(g2[i], fp);
//         element_init(g3[i], fp);
//     }
//
//     //---------------------------------
//     //  xi = xi^{(p-1)/6}
//     //---------------------------------
//     mpz_init(t);
//     mpz_init(u);
//     mpz_sub_ui(t, fp->order, 1);
//     mpz_fdiv_q_ui(t, t, 6);
//
//     element_init(xi, fp2);
//     element_init(tmp, fp2);
//     element_set_str(xi, "0 1");
//     element_pow(tmp, xi, t);
//
//     //---------------------------------
//     //  set gamma 1, 2, 3
//     //---------------------------------
//     element_set(g1[0], ((Element*)tmp->data)[0]);
//
//     for (i = 1; i < 5; i++) {
//         element_mul(g1[i], g1[i - 1], g1[0]);
//     }
//     for (i = 0; i < 5; i++) {
//         element_mul(g2[i], g1[i], g1[i]);
//     }
//     for (i = 0; i < 5; i++) {
//         element_mul(g3[i], g1[i], g2[i]);
//     }
//
//     pf->gamma1 = g1;
//     pf->gamma2 = g2;
//     pf->gamma3 = g3;
//
//     pf->glen1 = pf->glen2 = pf->glen3 = 5;
//
//     mpz_clear(t);
//     mpz_clear(u);
//
//     element_clear(xi);
//     element_clear(tmp);
// }
//
// //---------------------------------------------------------
// // precomputation for Fp12 operation
// //---------------------------------------------------------
// void bls509_fp24_precomp(Field f)
// {
//     field_precomp_p precomp = NULL;
//
//     precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));
//
//     precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
//
//     bls509_fp2_precomp_sqrt(precomp->ps, f);
//
//     precomp->pf = (field_precomp_frob_p)malloc(sizeof(struct ec_field_precomp_frob_st));
//
//     if (strcmp(f->field_name, "bls509_fp24a") == 0)
//     {
//         bls509_fp24_precomp_frob_beuchat(precomp->pf, f);
//     }
//
//     if (strcmp(f->field_name, "bls509_fp24b") == 0)
//     {
//         bls509_fp24_precomp_frob_aranha(precomp->pf, f);
//     }
//
//     f->precomp = (void *)precomp;
// }
//
// //---------------------------------------------------------
// // precomputation for Fp12 operation for pairing_init
// //---------------------------------------------------------
// void bls509_fp24_precomp_for_pairing_init(Field f)
// {
//     field_precomp_p precomp = NULL;
//
//     precomp = (field_precomp_p)malloc(sizeof(struct ec_field_precomp_st));
//
//     precomp->ps = (field_precomp_sqrt_p)malloc(sizeof(struct ec_field_precomp_sqrt_st));
//
//     bls509_fp2_precomp_sqrt_for_fp12init(precomp->ps, f);
//
//     precomp->pf = (field_precomp_frob_p)malloc(sizeof(struct ec_field_precomp_frob_st));
//
//     if (strcmp(f->field_name, "bls509_fp24a") == 0)
//     {
//         bls509_fp24_precomp_frob_beuchat(precomp->pf, f);
//     }
//
//     if (strcmp(f->field_name, "bls509_fp24b") == 0)
//     {
//         bls509_fp24_precomp_frob_aranha(precomp->pf, f);
//     }
//
//     f->precomp = (void *)precomp;
// }

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bls509_fp24_is_zero(const Element x)
{
    if (bls509_fp12_is_zero(rep1(x)))
    {
        if (bls509_fp12_is_zero(rep0(x))) {
            return TRUE;
        }
    }
    return FALSE;
}

int bls509_fp24_is_one(const Element x)
{
    if (bls509_fp12_is_zero(rep1(x)))
    {
        if (bls509_fp12_is_one(rep0(x))) {
            return TRUE;
        }
    }
    return FALSE;
}

int bls509_fp24_cmp(const Element x, const Element y)
{
    if (bls509_fp12_cmp(rep1(x), rep1(y)) == 0)
    {
        if (bls509_fp12_cmp(rep0(x), rep0(y)) == 0) {
            return 0;
        }
    }
    return 1;
}

int bls509_fp24_is_sqr(const Element x)
{
    Element *t = field(x)->base->tmp;

    if (element_is_zero(x)) {
        return FALSE;
    }

    bls509_fp12_sqr(t[0], rep0(x)); // t0 = x0^2
    bls509_fp12_sqr(t[1], rep1(x)); // t1 = x1^2

    bls509_fp12_gamma_mul(t[2], t[1]); // t1 = x1^2*xi
    bls509_fp12_add(t[0], t[0], t[2]);

    return bls509_fp12_is_sqr(t[0]);
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bls509_fp24_random(Element z)
{
    bls509_fp12_random(rep0(z));
    bls509_fp12_random(rep1(z));
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bls509_fp24_to_mpz(mpz_t a, const Element x)
{
    mpz_mul(a, rep(rep1(rep1(rep2(rep1(x))))), field(x)->base->base->base->base->order);   // a = rep121*p
    mpz_add(a, a, rep(rep0(rep1(rep2(rep1(x))))));   // a = a + rep111
    mpz_mul(a, a, field(x)->base->base->base->base->order);   //a = a*p
    mpz_add(a, a, rep(rep1(rep0(rep2(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep2(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep1(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep1(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep1(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep1(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep0(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep0(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep0(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep0(rep1(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep2(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep2(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep2(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep2(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep1(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep1(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep1(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep1(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep1(rep0(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep1(rep0(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep1(rep0(rep0(rep0(x))))));
    mpz_mul(a, a, field(x)->base->base->base->base->order);
    mpz_add(a, a, rep(rep0(rep0(rep0(rep0(x))))));
}

void bls509_fp24_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t s0;

    unsigned char b0[1536];
    mpz_t z;

    mpz_init(z);

    bls509_fp24_to_mpz(z, x);
    mpz_export(b0, &s0, 1, sizeof(*b0), 1, 0, z);

    memset(os, 0x00, 1536);

    memcpy(&os[1536 - (int)s0], b0, s0);

    (*size) = 1536;

    mpz_clear(z);
}

void bls509_fp24_from_oct(Element x, const unsigned char *os, const size_t size)
{
    mpz_t quo, rem;

    if (size < 1536) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_init(quo);
    mpz_init(rem);

    mpz_import(quo, size, 1, sizeof(*os), 1, 0, os);

    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep0(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep0(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep0(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep0(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep1(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep1(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep1(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep1(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep2(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep2(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep2(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep2(rep0(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep0(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep0(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep0(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep0(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep1(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep1(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep1(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep1(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep0(rep2(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep0(rep2(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep0(rep1(rep2(rep1(x))))), rem);
    mpz_tdiv_qr(quo, rem, quo, field(x)->base->base->base->base->order);
    mpz_set(rep(rep1(rep1(rep2(rep1(x))))), rem);

    mpz_clear(quo);
    mpz_clear(rem);
}
