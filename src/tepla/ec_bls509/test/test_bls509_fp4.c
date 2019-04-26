#include <assert.h>

#include "rdtsc.h"

#include "../ec_bls509_lcl.h"

#define N 1000
#define M 100

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
    Element d0, d1;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d0, f->base);
    element_init(d1, f->base);

    //--------------------
    //  add
    //--------------------
    element_set_str(a, "11CCC4D694F9682FBC508EB5410B250F61FE82E3456A0741AB5A06FA0BE6502C05193394437FCAA92D4E34E0039BCA9C0D750A890C6F7DF355ACFD5D7566CAE5 2BD3B6BAA573A9ABD63E569CC4785C9C27364E8BDA9F5BDE2CBF73737D95E4E8530ADA38E5F101D0F688E00C4451B033C44E25731587F4B5D2A37439603B788 F796244101D3FB235E368A61F3F446E22CEE87CB4FDAFFA1B00D99B52E87BEE947EA8737AC08571F416A8D2A3F3C90E52F1C91753DE53D7FE43C837E9B1A26E B4CB083CB25C894C40E615C84A766FC9B95D7674BCAE876022060512066D0D922183CC930430AD7AF079894EE578AD51508FBE8CBE0647CD94DE01FEF9F0B12");
    element_set_str(b, "F8288265535880757B43D7CD9621A1223EA4FC15403457624FED50758838755C9CE3DB97B78E0637F85AC0C245CA4E52DE4EEF68CCF9363B05E405B97D036C5 978BF42CCB802E37B229BF4C7F2B2CAB1A4EBFC453262C439BB4C0F27BD235CFFCA2E9253DF8DE013C44112BFD14567311B5BD4771B41FE96A4C184B4ACB207 94CE1806CEF80A505DD1299B53F51D85AC87203AA5E813829039A7F636D9F40CB21014EBB3B93EBDC4911353403C44D1819BDE22AA3E5FA0505BF51408D1B29 591CB55022A815B6CAA9FD6D4FA492837F00920E6797A848E5E126833C37DDB7BF9E07FE2E078C6EE056BD666BB9605B008CD19EC00C6CEB078EDF0428B1019");
    element_set_str(d, "BF9F5FCEAF5259B1716D97D20AC522ABA5E0A5B3C54CE2B9D6EDAFD845946C66C78EB8E42E00A1BDD08848BB63D9BAE4CF63C782AB17E5664CE2C2B54772EFF C35FAAE770F3D7E3886815E943A3894741850E502DC58821C8743465F9681AB84FADC35E23E9DFD232CCF138416606A6D603E2BA873C149F3CEF8C84AB0698F 370ECC47DD2F5BB3ED2888ADABDA94FB20C92370243B2A6111A7316D6458A73FD312402B9E3786D00945DA7663CB9887CA7C9F20FF4A6D1620C75FB717EEAEC 10DE7BD8CD5049F030B9013359A1B024D385E088324462FA907E72B9542A4EB49E121D491323839E9D0D046B551320DAC511C902B7E12B4B89C6CE10322A1B2B");

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

    // element_set_str(d, "EFDE90A01CC64464D95380934805A6F68C91A10653FD4EC75A4A7AA7B9F3B475B35C95BF22F7843414ABF787E9E6C797768F1D9F3258E59B1FABD24A6A06AF0 A86A6BB57788589EC83F782ECBCE05B306B49629D3A2988A534CF595E992FCA436DF4820FDDF83E5B629FB7AEA09BDB5EAF57F1F45362C052045689BB7B8C93 40DC67DB20E7975F8726C22AEF8B9D6F92E4A181492B51EBCF0D27924EA5B2E8E933D9C4D552A0B8ECAE647260423FFE9DD84852BBD0F67956D5CA75CAFDD39 1443F23958BCB32CBD2DB7DDFD4DC00901FD6103C76FC8AF00263AD29B028A457C0BD21CDCC81EF4B7A50B70E3F8A9BA34566417EE8F3608C59D54A4C34EDE37");

    element_set_str(d0, "EFDE90A01CC64464D95380934805A6F68C91A10653FD4EC75A4A7AA7B9F3B475B35C95BF22F7843414ABF787E9E6C797768F1D9F3258E59B1FABD24A6A06AF0 40DC67DB20E7975F8726C22AEF8B9D6F92E4A181492B51EBCF0D27924EA5B2E8E933D9C4D552A0B8ECAE647260423FFE9DD84852BBD0F67956D5CA75CAFDD39");
    element_set_str(d1, "A86A6BB57788589EC83F782ECBCE05B306B49629D3A2988A534CF595E992FCA436DF4820FDDF83E5B629FB7AEA09BDB5EAF57F1F45362C052045689BB7B8C93 1443F23958BCB32CBD2DB7DDFD4DC00901FD6103C76FC8AF00263AD29B028A457C0BD21CDCC81EF4B7A50B70E3F8A9BA34566417EE8F3608C59D54A4C34EDE37");

    element_set(((Element *)d->data)[0], d0);
    element_set(((Element *)d->data)[1], d1);

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_mul(c, a, b);
    }
    t2 = rdtsc();

    printf("element mul: %.2lf [clock]\n", (double)(t2 - t1) / N);

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
        element_inv(c, a);
    }
    t2 = rdtsc();

    printf("element inv: %.2lf [clock]\n", (double)(t2 - t1) / N);

    //--------------------
    //  pow
    //--------------------
    mpz_init(exp);

    mpz_set_str(exp, "AA4A2EE5234E0E95E8BB01F6B4A67F0EE8F2ADC1AA153C48D163AA85F3F534C", 16);

    element_pow(c, a, exp);
    element_set_str(d, "5A9AF4D04AEA90A6B3094478E0364E9DD9F0AE7919A3A4BC1E0C2562224C2D1A567B4B262BA3A49B7236CC5565DEDAED6BA7D26364BE1913099547E9A5410B0 440E6C9D20A9891EFDDBDB56F099A007B15E227D5A51ED333BB2641001517E50569CE98455293A13B60D0067AFEDF1C7C7DE215CACA59C231052B4F1B432B6A F47C3D9F8C367949BDEBC9C415D0EA63F5D1D32ACFCA387E496D49E7295C22AE7E6D497C7CBA4FD1A32CB016F2B90B038EAD76CC3B9AE4CF22249E5681B7C10 D11FA34995364937CC0A47495990295C0C12B0793F1489CE8D8B1DC3A079567D161207F13630E976F14D429C163FA81B6D40AB5BA945340AE82605B4C5F8BDE");

    assert(element_cmp(c, d) == 0);

    mpz_set(exp, f->order);

    for (i = 0; i < 50; i++)
    {
        element_random(a);

        element_pow(b, a, exp);

        assert(element_cmp(b, a) == 0);
    }

    t1 = rdtsc();
    for (i = 0; i < M; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with order: %.2lf [clock]\n", (double)(t2 - t1) / M);

    mpz_clear(exp);

    //--------------------
    //  clear
    //--------------------
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);

    element_clear(d0);
    element_clear(d1);
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

    for (i = 0; i < 50; i++)
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
    for (i = 0; i < M; i++) {
        element_sqrt(c, b);
    }
    t2 = rdtsc();
    printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / M);

    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(d);
}

//============================================
//   i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;
    char a_str[520];

    size_t blen;
    unsigned char b_str[256];

    Element a, b, c;

    element_init(a, f);
    element_init(b, f);
    element_init(c, f);

    for (i = 0; i < 100; i++)
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

    for (i = 0; i < 100; i++)
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

    // test for methods
    field_init(f, "bls509_fp4");
    test_feature(f);
    test_arithmetic_operation(f);
    // test_sqrt(f);
    test_io(f);

    field_clear(f);

    fprintf(stderr, "ok\n");

    return 0;
}
