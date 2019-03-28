#include <tepla/ec.h>

//----------------------------------------------
//  finite field : init clear
//----------------------------------------------
void ec_bls509_fp_new(Field f);
void ec_bls509_fp2_new(Field f);
void ec_bls509_fp4_new(Field f);
void ec_bls509_fp12_new(Field f);
void ec_bls509_fp24_new(Field f);

void ec_bls509_field_clear(Field f);

//----------------------------------------------
//  ellptic curve : init clear
//----------------------------------------------
void ec_bls509_fp_group_new(EC_GROUP ec);
void ec_bls509_tw_group_new(EC_GROUP ec);

void ec_bls509_group_clear(EC_GROUP ec);

//----------------------------------------------
//  pairing : init clear
//----------------------------------------------
void ec_bls509_pairing_new(EC_PAIRING p);

void ec_bls509_pairing_clear(EC_PAIRING p);
