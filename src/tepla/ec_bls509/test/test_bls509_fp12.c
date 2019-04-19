#include <assert.h>

#include "rdtsc.h"

#include "../ec_bls509_lcl.h"

#define N 100
#define M 10

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
    Element d1, d2, d3;
    Element f1, f2;
    Element g1, g2;
    Element h1, h2;

    mpz_t exp;

    //--------------------
    //  init
    //--------------------
    element_init(a, f);
    element_init(b, f);
    element_init(c, f);
    element_init(d, f);

    element_init(d1, f->base);
    element_init(d2, f->base);
    element_init(d3, f->base);

    element_init(f1, f->base->base);
    element_init(f2, f->base->base);

    element_init(g1, f->base->base);
    element_init(g2, f->base->base);

    element_init(h1, f->base->base);
    element_init(h2, f->base->base);


    //--------------------
    //  add
    //--------------------
    element_set_str(a, "21F369CDB0891084235F3C5C765072736D171E80AB47141AD34DAB4553FC099A005268582C46D9333193F863F440BA7E4E9739513FFEA76920835BCBAC4B992 348BDADC168B6B8BA7FD5C390C5B2966C241B4978A265C26F921C8C81540608085402EA34C093A2B749B2018388A6BBAA1E534550F476A087B281821066E205 152B898FC4D9A5D882169E2EB369F0CACBAEB17255C4140D731C92C62B24E786C2A9C0B448B801EFE1FBC21FEF75F73963D6C2E3D628EE22BE73D5AAF24249A2 4915D472D751D797D8F60E3895AB9C9E1D614DA72B0DB4E0E9A7F1AAC1FD2B78F276600A0C2F2553CB8D3C89D59E71D7B692A79CF36C2CAE802F441500ADEDE 10F68FD20F204ECE4E4BD94B85E4A590B615012AC25B41868F828E7B26E04E5E4F67EA5683946FA485B297FF87D5F35EFC9707A4E05E592D1EC0F60738926D87 F7F4EEFD651C0C6078317B1ACF2B51DC4343442F2D8B26FE4CB0BC681F9AB4B890946FF597874CD9013349D75A041B832CF40797DCADCDF8093BC9BA53AF7AD 106D25235AD83F48491EB1C21800453DC21EA19B833271B88E72E66D99972F5407E5009CF8CC5AD0DFF46D8BFA26584770BCE7A7BB68E45CF3EF1FA754D68F0B 14A25E5909BFE678C59C85B2357C26BF1B9C125EE555637BF957CEDE679F078D660C5704E01211B014CE965F8BCB31114ECC8EFB663547B08E9DEAAB7E75DAEA F6F2FB04B3D1905DDE4AE05F3C3A74F8A8E812CBD08DD331088B1BB505728E1E94604C4FF1BCED5317516312D1A5B57E0AB71045C6B0235A83D4F427BCCE159 CF59248032F02FE6E380947CB560BF2869942AAE95EE95E70296A1CD42DC7FED1FAB3A9D56A3CC9B47E44C03CD98C1D2375F3161FD7113A4C910C1A3E6358CC 58D62C79D699CD021E5D3ED3971F8C043385AD9B7F02CE7F7554D651148104E88454ACC9685177BDDF5D0DCCDC5356F6EC9C31F3BC31FE6A9FBC4AFE92914DD 5DC95E813999A7953BC3361D28AE751330847E961EA71F1E777ACCB120D36B47DC7CFE5FAFC2B1DDB04EDF3334B50392575775F1D4E995DAEDE182C4436D764");
    element_set_str(b, "11D0E8305F744165347E4B658D803ADAA28B1C528EA70015EC4AB28A66DE963BEEECAFE81B43317816C35D5CB2A8F7BEBF2F5E105F32B3B47007BDC10D8203E3 476FF5442BE0DDFE30B3636A8DC001461AF5AA450ECE4ECC5867E0BB17F65C5614497B9359F27BBCD8435CAA7B11836CF725313D12980FC1FF8CE2B4790E63E 11089178CDCC76A76C1EAFF0EBAA13DF01E03C9C22A80A873ACC3D02747510345E51C2701E1878B9416FE7E98A75CD67706F8455FF9626120D03263AACAF2FCC 49FB6C6BDBC7BF5D0A039E1301C4D61AB2FA446D100592E450E42247CA0E26AEEB03459B32AB2EC711D6839187F9106901BFF5E8F6710793EF7EA532E45DC9E 3562F5AC9D8B724A6F6E5F96DF9CDB54F1DE98D0768B758C03BABCEF38384B3EE061DEE0FB9D1DAE4BAF0D38BF82A343BC40DC6DBC3D06793BD6EC15175FEC0 A48B2B0BAB9E7DA958D062DC536C9F08158CBB394691C8679B25B312E36D74CD646962C91E691E9B036205CE8684DC8E36EADA0D8EF01197AF1750C50941232 1328745969D14DE37FDE1EA7BF53AF667DB8833A81B5BFBE9FCD68F49C25FE64432B5D084FA57C955AADE665D9C9B6FEF631A227192AEF8B0D11481A4CDF5170 87232B10EBA8AD6DBC48B84731F8BF9BB3E1E5407220DB6AC00581AB3436999D08AA6F7472E916DB5D7EA1313522C148939105C503E2A9507B5851250929208 15348A766B0C14C92A3B09B597323A3E7EEFBA2556F5F3F2DD380E09E64DF1854DDB61E0D0F7F94B0054060A9B9F2220E3F550BFF4D71982E18F75C88EBDD6DB 14A8D5075B345554D0C833FBED0923EC037AE9B308DA118F0E7E46391AB6C3B5DA39120B37124E209DC1202BB63B3509C759CA7C0F75BF20E6C8CB281FD4F6FC 13370C883F9158A3EF6EB4BE4924F8B8A4B5C63E70060513231CB62FBA9555AEE44B9AFE2B8FBDE2F7B88DE2D970D45AD6A54360051E710C731D4CD9E1AE94EF 882EC33013E1195B4B5227B0539CD281564754237F12F1DBF0500F6E60BCA3B694255FA527D53A26ED3C63B351FD432F80B6D4F90DEC8DA81F8162912B1842B");
    element_set_str(d, "13F01ECD3A7CD26D76B43F2B54E54201D95C8E3A995B7157997F8D3EBC1E56D58EF1D66D9E079F0B49DC9CE2F1ED0366A418D1A573329E2B020FF37DC846BD75 7BFBD020426C4989D8B0BFA39A1B2AACDD375EDC98F4AAF35189A9832D36BCD69989AA36A5FBB5E84CDE7CC2B39BEF27990A659221DF79CA7AB4FAD57F7C843 10DEC408936C51E3F1475B6AA55317B3020425C51B53A0087AFECEC4BF8966FFBE8CFD64EAB7D9B853A04DA90830F0CDE5E28A32673181342A39EA57E631A6C3 931140DEB31996F4E2F9AC4B977072B8D05B92143B1347C53A8C13F28C0B5227DD79A5A53EDA541ADD63C01B5D978240B8529D85E9DD34426FADE947E50BB7C 144CBF2CD8F905F2F542BF44F3DE73460532EAB7C9C3F8DF4FBE3A4A1A63D3123D6E0844934E417F6A6D88D313CE1D93385B156BBC222994B27E64C88A086C47 472AAA091D1DE04A0222B2A786892177A0237AD2A29506A2B9365F3D01FF1DCFCE1576C6F4665C6707DF899EC4DBBAE27DA3112E82C4AF85A48201A3D0F3734 E40427CC56FC28FCC0EDDB4DD9307AD744C5C8CA7CFB2EAFB564E5E55AC9CFCE8A1D7E5CC5936756AD6F79162353B73788ACCC7660640E75FC35633E8F60DD0 7BF3A0A1940A6B3A4731E81AEDAC5C20B4F68698F5EF2A6726E25F53AD1E06BD428783CAB28022CFADB24122D628952E9A1E25047E5DF44F5165E3016489A47 F4E6326B70F63330B31C5069134F4973DF37308B6E65299BAD6BEC1569489ABD4B2E0E653FB272F61FDBFDB56FEA9A5D63D04BCE2B488B7E88FB37D51CAE589 C49104F5F298DB742124A8EBE9E42E7BE89641495207C614BBDAF520ED3FAF949C53FF59063E9F98274088B8159ED53FC6C008AC0BF3D5A921CC5B4A5787D1D 36F184FDDC12AD8146695F688D604821C6358CECADDB36EE7880290EBCCD5420A22600B45FC346E05E3025F357B35F7570B4977D253FDF27BDBFFFC1217D721 E5F821B14D7AC0F087155DCD7C4B479486CBD2B99DBA10FA67CADC1F81900EFE70A25E04D797EC049D8B42E686B246C1D80E4AEAE2D623830D62E5556E85B8F");

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

    element_set_str(f1, "AF7F977BC13AB2AFD9454C10CE73F08B5321D0A64AF9E9084C3338CDC41E94FF7DEE750E1206842B671589B1F8981D44C240A634FADB418C4EF8AC41892C59F D5905636C0744BB8A1CA880947BE9868F4116FD46553B72B02E6FE355C5313A3D31F7C98B3C499411F395A17A2DB9CF2783B049A7BFFDCABFEDD5105F28CDC1");
    element_set_str(f2, "C3EDA2DAEF166A5C282F895056502E69D175A5D3BBC90C5C3DD76014B4BCBB60D405E2AAD4AF5FD8CD9A67A55DE2B31E6B7C30C1A72671AD4FABEF0D2CB592A 12724AA3487B9ECC3B72C86B0E6F518BBE6A3AD93D2783A08CCB0537835520A35E35C9274EBFBF4F60DD52063D366EC8C56056267D821E8883E90E105DF2B2EE");

    element_set_str(g1, "1008EF62EC4551BDA4067734B439F045C1485938F9D3BA8AFDF023D0C1A1713101F0BB7865D47281BFB5DAA354DBA843B20B2433A4B984F4A8E1BF36F0B006A8 1861623B85E134C6FF8C0651E6F5AACD15C6EC8C714DCC8E9CA762839BEFAED8D0F77FC011BDC205ACF2D8F775D9B1759C059F3171AE110CFF612F6CD275790");
    element_set_str(g2, "98216844E4AB109659F59286110D40004D879FD6B222171EA0CB2C19763E21D7B4D92A35E2CFD4E1A4A0E9C476F897C12AEBE1903804EC34F1CEB1B6305A75D D822AD83B831FE519FC265121F0D4DF5E2594E88CA124866DADD2E9A29F44B4B280E2F52B561C0D4F22656D973D193EFD9C59D7FC9B17D00F41EDBFAADB7A7B");

    element_set_str(h1, "5726ED4051F56DD868DDC322A522147A168874B80B850488521651EED72FCFEA69A0DB2EED49A99B338941ED5534790C72B86E6854A9517FC2B20C3B4A947B8 4BC588180007E051E968023FCAF977A8062E8390C762B46ABC74A6AC67BC3B89C6F2B6B867D39E9B5F5DA30876D93C1E11386677BB3E7E76B6E24FD9753DB67");
    element_set_str(h2, "FB690658DCFE5E9F0857DF742119F433B381BF4788FD95621A5DC8BCF261231A46D3C8EAFFFC91141CD5018AAA8051E1110EB7D48D326EEDEA8B117D0DDA84A 837E5919E6B78AA718E889D1519F2FC05F495147961E850305BB9E98638AAAA0ED8F0EFAD7A9C8DD90681DBAE6063BB5B54EC9E2236861602CE05B7E9738F13");

    element_set(((Element *)d1->data)[0], f1);
    element_set(((Element *)d1->data)[1], f2);

    element_set(((Element *)d2->data)[0], g1);
    element_set(((Element *)d2->data)[1], g2);

    element_set(((Element *)d3->data)[0], h1);
    element_set(((Element *)d3->data)[1], h2);

    element_set(((Element *)d->data)[0], d1);
    element_set(((Element *)d->data)[1], d2);
    element_set(((Element *)d->data)[2], d3);

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

    // char s[1560];
    // fprintf(stderr, "\n\n");
    // element_get_str(s, c);
    // fprintf(stderr, s);
    // fprintf(stderr, "\n\n");
    // element_get_str(s, d);
    // fprintf(stderr, s);
    // fprintf(stderr, "\n\n");

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

    mpz_set_str(exp, "103DC9103A8BDD2D08BA8E58EC4AD8AA108177A04CFBBAE6FB7CE07A88E298C3", 16);

    element_pow(c, a, exp);
    element_set_str(d, "11DB774CD07CCDD50F74AC8E8116611B0FC621231CE12151572B71D64065A348E6AF6CAFD988A40E3787A689F3D0DDB0F0905EE18A6DBCA743D2D3C2E067E076 B06EF3EF21E639D230FFE22B3AA5341CBFF2FCB87EA4F49837395B9C75507F94CA0806ED4067361D5568D15D941AF7246108429A62AC91AA49969878FBF2272 76DB41A127974B5B7F195844B608B07922487DED4D9825F31FFB7832C50CA5272E0E8F01D7C4F7984DE9D0EE1812F828FAE28A26BB7CB6172CEADBB851DA2FB 99B4CA1E1F7A7BDCCAE76490E9470AA1ADAFA78D633B50CE856EE68FF88C6B062DD54EEA776335831675C8075782E91CE3D3A58D1E90B696C7E2EF24CE59D43 79F015876A6036CDBAD503822292F32FABA8A2224AD6D9D7BE4F2AF4BF014195ACEF351B53200B594BA4EF61A6E822600CF0CCD4E32D25C2E5DB3E31CD7CF95 3B7FFB689A143E5DA542109E3F73001B1D7A10DA41589655D103B586091D6F05D7C7329E72922DB855B498523BEB9BEE38CBF3C880DDF52B0F9C87D1F701891 9CC35CA9CC03B3342C500EFCCFEDD2415EF9CB8B5D9768C4603D933E9B944373822A3A35FDFC5A8F96D10563A79FB430F1874EBFCE62FF4C560DD45507AD1F4 33DD5ECDAAFD3A9751AFC10CC7775DF14881A9F377F9388E687722D55C03B2AE7441AF21E5EADFF7362B62AE03B22F0046900ACD0EFD9C6A89A90ED740C7AC6 B00961EFAC63DF22FE5AC08A34ABB410A85CA013193138C97A577E14826807CC0829C92B43737F17E96414BE120C6910C1282F7B0152DA90FC911E0800F1CEA 90400FE463A9F130A2270DF7B74B602DB5B721A51EA7306179CE8FF9B5782099C8346BF507B0ED17DDEBF5ED62F9DFF3B51A81B35F017B558DF9D8FEC3F3EC3 34CDE45054DD39BF2AC4B014A33CAB8408380AFE5B550F06B126F90750846FA92AF3E02345DDBB4304E5569D30D3DB819AEACFEAFC2B0FBEFCA68BABA6A98B 125C22D72A47D5D9CA5132E34186FBD7ADB5AC87EEF8EE002259A9F1E4C9C67FFF9CB20808060972A2A24C0B8868691710B7E153A797252714B7120E1A32931D");

    assert(element_cmp(c, d) == 0);

    t1 = rdtsc();
    for (i = 0; i < N; i++) {
        element_pow(b, a, exp);
    }
    t2 = rdtsc();

    printf("element pow with torsion: %.2lf [clock]\n", (double)(t2 - t1) / N);

    mpz_set(exp, f->order);

    for (i = 0; i < 10; i++)
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

    element_clear(d1);
    element_clear(d2);
    element_clear(d3);

    element_clear(f1);
    element_clear(f2);
    element_clear(g1);
    element_clear(g2);
    element_clear(h1);
    element_clear(h2);
}

// //============================================
// //   test for sqrt
// //============================================
// void test_sqrt(Field f)
// {
//     int i;
//     unsigned long long int t1, t2;
//     Element a, b, c, d;
//
//     element_init(a, f);
//     element_init(b, f);
//     element_init(c, f);
//     element_init(d, f);
//
//     for (i = 0; i < 10; i++)
//     {
//         element_random(a);
//         element_sqr(b, a);
//
//         assert(element_is_sqr(b));
//
//         element_sqrt(c, b);
//         element_sqr(d, c);
//
//         assert(element_cmp(d, b) == 0);
//     }
//
//     t1 = rdtsc();
//     for (i = 0; i < N; i++) {
//         element_is_sqr(b);
//     }
//     t2 = rdtsc();
//
//     printf("element is sqr: %.2lf [clock]\n", (double)(t2 - t1) / N);
//
//     t1 = rdtsc();
//     for (i = 0; i < M; i++) {
//         element_sqrt(c, b);
//     }
//     t2 = rdtsc();
//
//     printf("element sqrt: %.2lf [clock]\n", (double)(t2 - t1) / M);
//
//     element_clear(a);
//     element_clear(b);
//     element_clear(c);
//     element_clear(d);
// }
//
// //============================================
// //   Frobenius Map \phi_p
// //============================================
// void test_frob(Field f)
// {
//     int i;
//     unsigned long long int t1, t2;
//     mpz_t p, p2, p3;
//     Element a, b, c;
//
//     mpz_init(p);
//     mpz_init(p2);
//     mpz_init(p3);
//
//     mpz_set(p, *field_get_char(f));
//     mpz_mul(p2, p, p);
//     mpz_mul(p3, p2, p);
//
//     element_init(a, f);
//     element_init(b, f);
//     element_init(c, f);
//
//     for (i = 0; i < 100; i++)
//     {
//         element_random(a);
//         element_pow(b, a, p);
//         bls509_fp12_frob_p(c, a);
//
//         assert(element_cmp(b, c) == 0);
//     }
//
//     t1 = rdtsc();
//     for (i = 0; i < N; i++) {
//         bls509_fp12_frob_p(c, a);
//     }
//     t2 = rdtsc();
//
//     printf("element frob p: %.2lf [clock]\n", (double)(t2 - t1) / N);
//
//     for (i = 0; i < 100; i++)
//     {
//         element_random(a);
//         element_pow(b, a, p2);
//         bls509_fp12_frob_p2(c, a);
//
//         assert(element_cmp(b, c) == 0);
//     }
//
//     t1 = rdtsc();
//     for (i = 0; i < N; i++) {
//         bls509_fp12_frob_p2(c, a);
//     }
//     t2 = rdtsc();
//
//     printf("element frob p2: %.2lf [clock]\n", (double)(t2 - t1) / N);
//
//     for (i = 0; i < 100; i++)
//     {
//         element_random(a);
//         element_pow(b, a, p3);
//         bls509_fp12_frob_p3(c, a);
//
//         assert(element_cmp(b, c) == 0);
//     }
//
//     t1 = rdtsc();
//     for (i = 0; i < N; i++) {
//         bls509_fp12_frob_p3(c, a);
//     }
//     t2 = rdtsc();
//
//     printf("element frob p3: %.2lf [clock]\n", (double)(t2 - t1) / N);
//
//     mpz_clear(p);
//     mpz_clear(p2);
//     mpz_clear(p3);
//
//     element_clear(a);
//     element_clear(b);
//     element_clear(c);
// }

//============================================
//   i/o test
//============================================
void test_io(Field f)
{
    int i;
    unsigned long long int t1, t2;
    char a_str[1560];

    size_t blen;
    unsigned char b_str[768];

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
//  main program
//============================================
int main(void)
{
    Field f;

    // test for methods
    field_init(f, "bls509_fp12");
    test_feature(f);
    test_arithmetic_operation(f);
    // test_sqrt(f);
    // test_frob(f);
    test_io(f);

    field_clear(f);

    fprintf(stderr, "ok\n");

    return 0;
}
