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
    "36854167537133870167810883151830777579616207957825464098945783786886075923783763188360549476" \
    "76345821548104185464507 "                                                                     \
    "13395065449444764730204713799419212215849338759383496204265437364165114239563335064727246553" \
    "53366534992391756441569"
#define HGEN                                                                                       \
    "1 "                                                                                           \
    "35270106958746661818713911601106014489002995279277524021990864423979378573571502687334760034" \
    "3865175952761926303160 "                                                                      \
    "30591443442442137099712598147537816369864703254766475586593732062916353247689584324335095631" \
    "04347017837885763365758 "                                                                     \
    "19851506022872919355680545211771716383008689782156557308593786650663447263738237184238691042" \
    "63333984641494340347905 "                                                                     \
    "92755366549233245574720196577603788075774019345359297002502797879397687700267556498094928972" \
    "7957565575433344219582"
#define SEED1 "36927108280682610572751541835651759786319781688803963047055892652598699269454"
#define SEED2 "48642930395851658792769562567117246310440249719334092526800582777656115356016"
#define SEED3 "9337391987890516768459279655811256076141636705119217609748065801581636739148"
#endif

#include "CUnit/Basic.h"
#include "../src/common/sha256.c"
#include <gmp.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

#ifndef IS_MAC_OS
#include <sys/sysinfo.h>
#endif

int logN;
int Nb;
int Mc;

static int M;
int N;
int nPublic;
int nConst;
int setParams;

mclBnFr* uw;

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

#include "gro16.h"

void bulletproof_prove(unsigned char* si[]);
int bulletproof_verify();
void bulletproof_save();
void bulletproof_read();
static inline void bulletproof_init(int Nb_set, int Mc_set);
static inline void bulletproof_get_context(context* ctx);
static inline void bulletproof_user_gammas(int val);
void init_setup(void* circuit);
setup_keys perform_setup(void* circuit);
void init_prover(void* circuit, proving_key* pk);
proof generate_proof(void* circuit, proving_key* pk);
int verify_proof(void* circuit, proof* p, verifying_key* vk);

#include "../src/bulletproofs/bulletproofs.c"
#include "../src/zpie.c"