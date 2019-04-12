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
    element_set_str(a, "500AB993D30A0A9905DE403C439132E893498DB76A786FE53763FE1066E2559ED6AB2ADEE510EDE6FF2B7750EF6C72723C644BC507576C01EB74BF7CD52131F 13836C4BB937ED5BD0BBA5EBDEA605098C234F7D19A8C33BAC191FCD4CB61506FB106A60CD21FA5F9FD213E329EC2DC36BC130AC9962C547E0910F8CC1307165 13A25FBCB2426F572D176AF9D408334477DDC817605A1DBB7D3C39CC4CA03A026114C75D714EAEC2913FE15D45ADBD04B8C004E98321585D2FFF6B94F39B2465 473A1CB40908E80144E41BC75E4FF9293E90C1AA82C8256C39ED3D60AECD3B6D5B2099FBBE7380D828C231048EA64B7311A9EB9824548AE2705ED51AF61CE30");
    element_set_str(b, "904E39AE9B6B81592DE83DA933DB3E89BF93E5322545AF805B52C367ECA7D7EE4BA03C2B5828A6F0667D9D3AE1CC429E781BDE5E6DD742DCEE20D6529E1C6C7 8412D9214AA88C4C5E0BF5FF4E124813C4D4025266282A3431BDD339D2217CED3ACEA88ADE1CD4BC02D79776773D4361F6CC22042CF463C5E09753CC3D54107 151C1F95ED8A6F35A4513174361EC9D4CA8AB8A74AD8EE10286282B05A19FC050C8ECA4A405535C5723A77E19E783C6450B5B3BEF3722F6D2FC8EE410926037F 1435E448A7D7AA9AEC7E2187F3BA5C99A3ED94BEB75BDC0F89915BD15318D8663AFEAD4946F301978FE6BE2CA0C7B18E1A62A76B3B788DECF82CBC43A12BD2E1");
    element_set_str(d, "E058F3426E758BF233C67DE5776C717252DD72E98FBE1F6592B6C178538A2D8D224B670A3D3994D765A9148BD138B510B4802A23752EAEDED99595CF733D9E6 66F42DDCEA8AB8499AE7296D9C63C93FCE5C758E2F2C752BC4AFBFD09C79C1A6C4ECF29FEEB26BA903430FA1FA52E269CCA35C56DA478839D5D733BCC45DFC1 13692852A09313F0D47AA9B91066102276DDB8754E1A8D3F72B4BB78C6A9A54C0B350BE8358B439733AEFCDE726B25961B11FBA10805F4C9BE8B484844015539 3542F13E92E6E7F03DE708F6FDE6F356C4BD890026FDFDA1A462EA37DF51B61AE42312986C198B442A784DC77F742725D19891D4F30439A7DF5980797CDCE66");

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

    element_set_str(d, "21628B383269A698B7FADDD5D525C0AA747F10BD3DC795DC3D0E249A5E42524DE6495B778255D353BC7741297D25DCB39E991D6DD6F4590C3D417987C9C9E15 B5A2355FF505892D77882A3A555608486B776BA33527B029FCE7F33BA5BCC5C090D439265E596D5C5A01419F8F7E21DBB0331C4E3E67F3D64F4885FF9E97553 BC5F211DAAA3B1DAD6E4788427618B4206D955C024F3EAAB392656E3E10F187257227A6273663F378F78CAD5DA71333B3DEAAE95C1A44ED9E5DBFCC418B62 152016B7323C9A0E460D98BC4B63F3749E963D7CC5FC03B1D8727D8BC316E21CF21A67E41819812F151C1B59C28EDC1EBA9513F8A62090C0BCC0E1A8ABAE1BF6");
    // element_set_str(d0, "21628B383269A698B7FADDD5D525C0AA747F10BD3DC795DC3D0E249A5E42524DE6495B778255D353BC7741297D25DCB39E991D6DD6F4590C3D417987C9C9E15 BC5F211DAAA3B1DAD6E4788427618B4206D955C024F3EAAB392656E3E10F187257227A6273663F378F78CAD5DA71333B3DEAAE95C1A44ED9E5DBFCC418B62");
    // element_set_str(d1, "B5A2355FF505892D77882A3A555608486B776BA33527B029FCE7F33BA5BCC5C090D439265E596D5C5A01419F8F7E21DBB0331C4E3E67F3D64F4885FF9E97553 152016B7323C9A0E460D98BC4B63F3749E963D7CC5FC03B1D8727D8BC316E21CF21A67E41819812F151C1B59C28EDC1EBA9513F8A62090C0BCC0E1A8ABAE1BF6");
    //
    // element_set(((Element *)d->data)[0], d0);
    // element_set(((Element *)d->data)[1], d1);

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
    element_set_str(d, "A9942F0FA93B9133503E9940B6FEC891AE9851AF39A3616ABFE403B2EB5F1C0B10F4D07AB77F762EDDCD8449CCD9068C4DA8187FC92B7E5C0B2646D79C26BCB C5D7D4F3F7DF240F94B5A5BA7AED19B589A75685515EDDBF69E03C2137F6D984A94862BA6D1F5742AEF0F156E3F6F42E9D82D4FD00B1C940B06DC1FD64FCE50 F3A0377C8DAD6A78ED67A850E906E75645175CCA6CF02D3C963FA9F259DD43E8ABE61A4B2047C3E989D87A6B49575A6F685738CD1C239AAFDFFB6AC68185F4F C09915FDFC06185B1AB2A60E60A63298DB569B0B2C69F8CB6D0D1DE8B39E3E4EB88B33D9590232A4CE04E93A034767EBF24F5F7B9FD845EB20C4005CDC0E084");

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
    char a_str[390];

    size_t blen;
    unsigned char b_str[192];

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
    test_sqrt(f);
    test_io(f);

    field_clear(f);

    fprintf(stderr, "ok\n");

    return 0;
}
