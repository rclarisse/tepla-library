//==============================================================
//  prime field ( bls509_fp ) implementation with GMP
//--------------------------------------------------------------
//  2019.03.25 created by rclarisse
//==============================================================
#include <pthread.h>
#include "ec_bls509_lcl.h"

#define rep(x)  (*((mpz_t *)x->data))
#define order(x) (x->field->order)
#define field(x) (x->field)

//-------------------------------------------
//  initialization, clear, set
//-------------------------------------------
void bls509_fp_init(Element x)
{
    x->data = (void *)malloc(sizeof(mpz_t));

    if (x->data == NULL) {
        fprintf(stderr, "fail: allocate in fp init\n");
        exit(100);
    }

    mpz_init(rep(x));
}

void bls509_fp_clear(Element x)
{
    if (x->data != NULL)
    {
        mpz_clear(rep(x));
        free(x->data);
        x->data = NULL;
    }
}

//-------------------------------------------
//  set value
//-------------------------------------------
void bls509_fp_set(Element x, const Element y)
{
    mpz_set(rep(x), rep(y));
}

void bls509_fp_set_str(Element x, const char *s)
{
    mpz_set_str(rep(x), s, 16);
    mpz_mod(rep(x), rep(x), order(x));
}

void bls509_fp_get_str(char *s, const Element x)
{
    mpz_get_str(s, 16, rep(x));
}

void bls509_fp_set_zero(Element x)
{
    mpz_set_ui(rep(x), 0);
}

void bls509_fp_set_one(Element x)
{
    mpz_set_ui(rep(x), 1);
}

//-------------------------------------------
//  arithmetic operation
//-------------------------------------------
void bls509_fp_add(Element z, const Element x, const Element y)
{
    mpz_add(rep(z), rep(x), rep(y));
    if (mpz_cmp(rep(z), order(z)) >= 0)
    {
        mpz_sub(rep(z), rep(z), order(z));
    }
}

void bls509_fp_addn(Element z, const Element x, const Element y)
{
    mpz_add(rep(z), rep(x), rep(y));
}

void bls509_fp_addp(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), order(x));
}

void bls509_fp_add_one(Element z, const Element x)
{
    mpz_add_ui(rep(z), rep(x), 1);
}

void bls509_fp_neg(Element z, const Element x)
{
    mpz_sub(rep(z), order(z), rep(x));
}

void bls509_fp_sub(Element z, const Element x, const Element y)
{
    mpz_sub(rep(z), rep(x), rep(y));
    if (mpz_sgn(rep(z)) < 0)
    {
        mpz_add(rep(z), rep(z), order(z));
    }
}

void bls509_fp_subn(Element z, const Element x, const Element y)
{
    mpz_sub(rep(z), rep(x), rep(y));
}

void bls509_fp_mul(Element z, const Element x, const Element y)
{
    mpz_mul(rep(z), rep(x), rep(y));
    mpz_mod(rep(z), rep(z), order(z));
}

void bls509_fp_muln(Element z, const Element x, const Element y)
{
    mpz_mul(rep(z), rep(x), rep(y));
}

void bls509_fp_mulc(Element z, const Element x, const mpz_t c)
{
    mpz_mul(rep(z), rep(x), c);
}

void bls509_fp_div2(Element z, const Element x)
{
    mpz_set(rep(z), rep(x));
    if (mpz_odd_p(rep(z)))
    { mpz_add(rep(z), rep(z), order(z)); }
    mpz_fdiv_q_2exp(rep(z), rep(z), 1);
}

void bls509_fp_inv(Element z, const Element x)
{
    mpz_invert(rep(z), rep(x), order(z));
}

void bls509_fp_mod(Element z, const Element x)
{
    mpz_mod(rep(z), rep(x), order(z));
}

void bls509_fp_dob(Element z, const Element x)
{
    bls509_fp_add(z, x, x);
}

void bls509_fp_tri(Element z, const Element x)
{
    Element *t = field(z)->tmp;

    bls509_fp_add(t[0], x, x);
    bls509_fp_add(z, t[0], x);
}

void bls509_fp_sqr(Element z, const Element x)
{
    mpz_mul(rep(z), rep(x), rep(x));
    mpz_mod(rep(z), rep(z), order(z));
}

void bls509_fp_pow(Element z, const Element x, const mpz_t exp)
{
    mpz_powm(rep(z), rep(x), exp, order(z));
}

int bls509_fp_sqrt(Element z, const Element x)
{
    mpz_t exp;

    if (!bls509_fp_is_sqr(x)) {
        return FALSE;
    }

    mpz_init_set(exp, order(x));

    mpz_add_ui(exp, exp, 1);
    mpz_fdiv_q_2exp(exp, exp, 2);

    bls509_fp_pow(z, x, exp);

    mpz_clear(exp);

    return TRUE;
}

void bls509_fp_OP1_1(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP1_1);
}

void bls509_fp_OP1_2(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP1_2);
}

void bls509_fp_OP2(Element z, const Element x)
{
    mpz_add(rep(z), rep(x), field(x)->OP2);
}

//-------------------------------------------
//  comparison operation
//-------------------------------------------
int bls509_fp_is_zero(const Element x)
{
    return (mpz_sgn(rep(x)) == 0);
}

int bls509_fp_is_one(const Element x)
{
    return (mpz_cmp_ui(rep(x), 1) == 0);
}

int bls509_fp_is_sqr(const Element x)
{
    return (mpz_legendre(rep(x), order(x)) == 1);
}

int bls509_fp_cmp(const Element x, const Element y)
{
    return mpz_cmp(rep(x), rep(y));
}

//-------------------------------------------
//  general function for is sqr
//-------------------------------------------
int bls509_fp_is_sqr_general(const Element x)
{
    int hr = FALSE;

    mpz_t q;
    Element t;

    if (element_is_zero(x)) {
        return FALSE;
    }

    mpz_init(q);
    element_init(t, field(x));

    mpz_sub_ui(q, order(x), 1);
    mpz_tdiv_q_2exp(q, q, 1);

    element_pow(t, x, q);

    hr = element_is_one(t);

    mpz_clear(q);
    element_clear(t);

    return hr;
}

//-------------------------------------------
//  generate random element
//-------------------------------------------
void bls509_fp_random(Element z)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    static gmp_randstate_t *state = NULL;

    pthread_mutex_lock(&mutex);
    if (state == NULL)
    {
        state = (gmp_randstate_t *)malloc(sizeof(gmp_randstate_t));
        gmp_randinit_default(*state);
        gmp_randseed_ui(*state, (int)time(NULL));
    }
    mpz_urandomm(rep(z), *state, order(z));
    pthread_mutex_unlock(&mutex);
}

//-------------------------------------------
//  i/o operation (octet string)
//-------------------------------------------
void bls509_fp_to_oct(unsigned char *os, size_t *size, const Element x)
{
    size_t stmp;

    unsigned char ostmp[64];

    mpz_export(ostmp, &stmp, 1, sizeof(*ostmp), 1, 0, rep(x));

    memset(os, 0x00, 64);

    memcpy(&(os[64 - (int)stmp]), ostmp, stmp);

    (*size) = 64;
}

void bls509_fp_from_oct(Element x, const unsigned char *os, const size_t size)
{
    if (size < 64) {
        fprintf(stderr, "error: please set up the enought buffer for element\n");
        exit(300);
    }

    mpz_import(rep(x), size, 1, sizeof(*os), 1, 0, os);
}
