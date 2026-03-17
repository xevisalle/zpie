/* SPDX-License-Identifier: GPL-3.0
 *
 * zpie.h - Zero-Knowledge Proofs in Embedded Systems (public API)
 *
 * Author: xevisalle
 */

#ifndef ZPIE_H
#define ZPIE_H

#include <mcl/bn_c384_256.h>

/* ------------------------------------------------------------------------- */
/* Public types                                                              */
/* ------------------------------------------------------------------------- */

typedef struct
{
    int index;
} zpie_element;

typedef struct
{
    mclBnFr* uwProof;
    mclBnG1 piA, piC;
    mclBnG2 piB2;
} zpie_proof;

typedef struct
{
    mclBnFr* constants;
    mclBnGT alphabetaT;
    mclBnG2 gamma2;
    mclBnG2 delta2;
    mclBnG1* vk1;
} zpie_verifying_key;

typedef struct
{
    int Ne;
    mclBnFr* wM;

    int qap_size;
    int* LRO;
    mclBnFr* LRO_constants;

    mclBnG1 alpha1;
    mclBnG1 beta1;
    mclBnG2 beta2;
    mclBnG1 delta1;
    mclBnG2 delta2;

    mclBnG1* A1;
    mclBnG1* B1;
    mclBnG2* B2;
    mclBnG1* pk1;
    mclBnG1* xt1;
    mclBnG1* xt1_rand;
} zpie_proving_key;

typedef struct
{
    zpie_proving_key pk;
    zpie_verifying_key vk;
} zpie_setup_keys;

/* ------------------------------------------------------------------------- */
/* Groth'16 API                                                              */
/* ------------------------------------------------------------------------- */

void zpie_init(zpie_element* toAdd);
void zpie_init_array(zpie_element* toAdd, int size);
void zpie_init_public(zpie_element* toAdd);
void zpie_input(zpie_element* var, char* val);
void zpie_mul(zpie_element* oo, zpie_element* lo, zpie_element* ro);
void zpie_assert_equal(zpie_element* lo, zpie_element* ro);
void zpie_addmul(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* ro);
void zpie_add3mul(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* lo3, zpie_element* ro);
void zpie_addmuladd(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* ro1, zpie_element* ro2);
void zpie_add3muladd3(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* lo3, zpie_element* ro1,
                      zpie_element* ro2, zpie_element* ro3);
void zpie_addsmul(zpie_element* oo, int* size, zpie_element* los, zpie_element* ro);
void zpie_addmul_constants(zpie_element* oo, int* lc1, zpie_element* lo1, int* lc2, zpie_element* lo2, int* rc,
                           zpie_element* ro);
void zpie_mul_constants(zpie_element* oo, int* lc, zpie_element* lo, int* rc, zpie_element* ro);
void zpie_mul_big_constants(zpie_element* oo, mclBnFr* lc, zpie_element* lo, mclBnFr* rc, zpie_element* ro);

void zpie_perform_setup(zpie_setup_keys* keys, void* circuit);
void zpie_generate_proof(zpie_proof* p, void* circuit, zpie_proving_key* pk);
int zpie_verify_proof(void* circuit, zpie_proof* p, zpie_verifying_key* vk);

void zpie_store_setup(zpie_setup_keys* keys);
void zpie_read_setup(zpie_setup_keys* keys, void* circuit);
void zpie_store_proof(zpie_proof* p);
void zpie_read_proof(zpie_proof* p);

/* ------------------------------------------------------------------------- */
/* Bulletproofs API                                                          */
/* ------------------------------------------------------------------------- */

struct context;

void zpie_bulletproof_prove(unsigned char* si[]);
int zpie_bulletproof_verify();
void zpie_bulletproof_save();
void zpie_bulletproof_read();
static inline void zpie_bulletproof_init(int Nb_set, int Mc_set);
static inline void zpie_bulletproof_get_context(struct context* ctx);
static inline void zpie_bulletproof_user_gammas(int val);

#endif /* ZPIE_H */
