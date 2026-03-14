/* SPDX-License-Identifier: GPL-3.0
 *
 * zpie.h - Zero-Knowledge Proofs in Embedded Systems
 *
 * Author: xevisalle
 */

#include <mcl/bn_c384_256.h>

#ifdef BN128
#define USEDCURVE MCL_BN_SNARK1
#define BITS 254
#define SIZE_FR 32
#define SIZE_G1 32
#define SIZE_G2 64
#define SIZE_GT 384
#define PRIMESTR "21888242871839275222246405745257275088548364400416034343698204186575808495617"
#define GROUPGEN 5
#define GGEN "1 1 2"
#define HGEN                                                                                       \
    "1 10857046999023057135944570762232829481370756359578518086990519993285655852781 "             \
    "11559732032986387107991004021392285783925812861821192530917403151452391805634 "               \
    "8495653923123431417604973247489272438418190587263600148770280649306958101930 "                \
    "4082367875863433681332203403145435568316851327593401208105741076214120093531"
#define SEED1 "12240074246416256392170098003078167155411687643283363511613846943233267248740"
#define SEED2 "343068200766678522613540249871646276936501879543207507284832640577391066390"
#define SEED3 "2638070508314633938838656314678915510188754475627216070946991168461793987539"
#elif BLS12_381
#define USEDCURVE MCL_BLS12_381
#define BITS 255
#define SIZE_FR 32
#define SIZE_G1 48
#define SIZE_G2 96
#define SIZE_GT 576
#define PRIMESTR "52435875175126190479447740508185965837690552500527637822603658699938581184513"
#define GROUPGEN 7
#define GGEN                                                                                       \
    "1 "                                                                                           \
    "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054"     \
    "9476"                                                                                         \
    "76345821548104185464507 "                                                                     \
    "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724"     \
    "6553"                                                                                         \
    "53366534992391756441569"
#define HGEN                                                                                       \
    "1 "                                                                                           \
    "3527010695874666181871391160110601448900299527927752402199086442397937857357150268733476"     \
    "0034"                                                                                         \
    "3865175952761926303160 "                                                                      \
    "3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509"     \
    "5631"                                                                                         \
    "04347017837885763365758 "                                                                     \
    "1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869"     \
    "1042"                                                                                         \
    "63333984641494340347905 "                                                                     \
    "9275536654923324557472019657760378807577401934535929700250279787939768770026755649809492"     \
    "8972"                                                                                         \
    "7957565575433344219582"
#define SEED1 "36927108280682610572751541835651759786319781688803963047055892652598699269454"
#define SEED2 "48642930395851658792769562567117246310440249719334092526800582777656115356016"
#define SEED3 "9337391987890516768459279655811256076141636705119217609748065801581636739148"
#endif

#include <limits.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#ifndef IS_MAC_OS
#include <sys/sysinfo.h>
#endif

/* ------------------------------------------------------------------------- */
/* Protocol structs                                                          */
/* ------------------------------------------------------------------------- */

typedef struct
{
    mclBnFr* gammas;
    mclBnG1* V;
    mclBnG1 G;
    mclBnG1 H;
} context;

typedef struct
{
    mclBnFr* uwProof;
    mclBnG1 piA, piC;
    mclBnG2 piB2;
} proof;

typedef struct
{
    mclBnFr* constants;
    mclBnGT alphabetaT;
    mclBnG2 gamma2;
    mclBnG2 delta2;
    mclBnG1* vk1;
} verifying_key;

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
} proving_key;

typedef struct
{
    proving_key pk;
    verifying_key vk;
} setup_keys;

typedef struct
{
    int index;
} element;

struct Trapdoor
{
    mclBnFr alpha;
    mclBnFr beta;
    mclBnFr gamma;
    mclBnFr delta;
    mclBnFr x;
};

struct Sigma1
{
    mclBnG1 alpha;
    mclBnG1 beta;
    mclBnG1 delta;
    mclBnG1* A;
    mclBnG1* B;
    mclBnG1* vk;
    mclBnG1* pk;
    mclBnG1* xt;
};

struct Sigma2
{
    mclBnG2 beta;
    mclBnG2 gamma;
    mclBnG2 delta;
    mclBnG2* B;
};

struct mulExpResult
{
    mclBnG1 htdelta;
    mclBnG1 uwA1;
    mclBnG1 uwB1;
    mclBnG2 uwB2;
    mclBnG1 uwC1;
};

/* ------------------------------------------------------------------------- */
/* Helpers                                                                   */
/* ------------------------------------------------------------------------- */

extern int test_no_rand;
extern element one, oneNeg, c_mimc[91];
extern mclBnFr* uw;
extern int setParams;

extern int M;
extern int N;
extern int nPublic;
extern int nConst;
extern int setParams;

extern int bench;
extern int logs;
extern int prover;
extern int lro_const_total;
extern int lro_constants_n;
extern int constant_n;
extern int un;
extern int wn;

/* ------------------------------------------------------------------------- */
/* Groth'16 API                                                              */
/* ------------------------------------------------------------------------- */

void init(element* toAdd);
void init_array(element* toAdd, int size);
void init_public(element* toAdd);
void input(element* var, char* val);
void mul(element* oo, element* lo, element* ro);
void assert_equal(element* lo, element* ro);
void addmul(element* oo, element* lo1, element* lo2, element* ro);
void add3mul(element* oo, element* lo1, element* lo2, element* lo3, element* ro);
void addmuladd(element* oo, element* lo1, element* lo2, element* ro1, element* ro2);
void add3muladd3(element* oo, element* lo1, element* lo2, element* lo3, element* ro1, element* ro2,
                 element* ro3);
void addsmul(element* oo, int* size, element* los, element* ro);
void addmul_constants(element* oo, int* lc1, element* lo1, int* lc2, element* lo2, int* rc,
                      element* ro);
void mul_constants(element* oo, int* lc, element* lo, int* rc, element* ro);
void mul_big_constants(element* oo, mclBnFr* lc, element* lo, mclBnFr* rc, element* ro);

setup_keys perform_setup(void* circuit);
proof generate_proof(void* circuit, proving_key* pk);
int verify_proof(void* circuit, proof* p, verifying_key* vk);

void init_circuit(void* circuit);
void init_setup(void* circuit);
void init_prover(void* circuit, proving_key* pk);
void store_setup(setup_keys* keys);
proof read_proof();
setup_keys read_setup(void* circuit);
void store_proof(proof* p);

/* ------------------------------------------------------------------------- */
/* Bulletproofs API                                                          */
/* ------------------------------------------------------------------------- */

void bulletproof_prove(unsigned char* si[]);
int bulletproof_verify();
void bulletproof_save();
void bulletproof_read();
static inline void bulletproof_init(int Nb_set, int Mc_set);
static inline void bulletproof_get_context(context* ctx);
static inline void bulletproof_user_gammas(int val);
