#include <mcl/bn_c384_256.h>

#ifdef BN128
    #define USEDCURVE MCL_BN_SNARK1
    #define BITS 254
    #define PRIMESTR "21888242871839275222246405745257275088548364400416034343698204186575808495617"
    #define GROUPGEN 5
    #define GGEN "1 1 2"
    #define HGEN "1 10857046999023057135944570762232829481370756359578518086990519993285655852781 11559732032986387107991004021392285783925812861821192530917403151452391805634 8495653923123431417604973247489272438418190587263600148770280649306958101930 4082367875863433681332203403145435568316851327593401208105741076214120093531"
    #define SEED1 "12240074246416256392170098003078167155411687643283363511613846943233267248740"
    #define SEED2 "343068200766678522613540249871646276936501879543207507284832640577391066390"
    #define SEED3 "2638070508314633938838656314678915510188754475627216070946991168461793987539"
#elif BLS12_381
    #define USEDCURVE MCL_BLS12_381
    #define BITS 255
    #define PRIMESTR "52435875175126190479447740508185965837690552500527637822603658699938581184513"
    #define GROUPGEN 7
    #define GGEN "1 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569"
    #define HGEN "1 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582"
    #define SEED1 "36927108280682610572751541835651759786319781688803963047055892652598699269454" 
    #define SEED2 "48642930395851658792769562567117246310440249719334092526800582777656115356016" 
    #define SEED3 "9337391987890516768459279655811256076141636705119217609748065801581636739148"  
#endif

#include <stdio.h>
#include <gmp.h>
#include <sys/stat.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <math.h> 
#include "sha256.c" 

int logN;
int Nb;
int Mc; 

int M;
int N;
int nPublic = 0;
int setParams;

mpz_t *uw;

typedef struct Context
{
    mclBnFr *gammas;
    mclBnG1 *V;
    mclBnG1 G;
    mclBnG1 H;
} context;

#include "gro16.h"

void bulletproof_prove(unsigned char *si[]);
int bulletproof_verify();
void bulletproof_save();
void bulletproof_read();
static inline void bulletproof_init(int Nb_set, int Mc_set);
static inline void bulletproof_get_context(context *ctx);
static inline void bulletproof_user_gammas(int val);
void init_setup();
void perform_setup();
void init_prover();
void generate_proof();
void init_verifier();
int verify_proof();

#include "zpie.c"
#include "bulletproofs.c"