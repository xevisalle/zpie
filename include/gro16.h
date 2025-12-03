#ifndef GRO16_H
#define GRO16_H

int bench;

typedef struct
{
    int index;
} element;

element one, oneNeg, c_mimc[91];

int logs;
int test_no_rand;

void element_log(char* text, element* oo);
void init(element* toAdd);
void circuit();
void init_array(element* toAdd, int size);
void init_public(element* toAdd);
void setPublic(element* set);
void input(element* var, char* val);
void mul(element* oo, element* lo, element* ro);
void assert_equal(element* lo, element* ro);
void addmul(element* oo, element* lo1, element* lo2, element* ro);
void submul(element* oo, element* lo1, element* lo2, element* ro);
void add3mul(element* oo, element* lo1, element* lo2, element* lo3, element* ro);
void addmuladd(element* oo, element* lo1, element* lo2, element* ro1, element* ro2);
void mpz_to_fr(mclBnFr* frFactor, mpz_t* convert);
void fr_to_mpz(mpz_t* convert, mclBnFr* frFactor);
int fr_cmp(mclBnFr* frFactor1, mclBnFr* frFactor2);
void sort_list(mpz_t* exp[], int heapsize);
void binarymaxheap(mpz_t* exp[], int i, int heapsize);

static mpz_t pPrime;
static gmp_randstate_t state;

int** L;
int** R;
int** O;

mpz_t* LRO_constants;

static mclBnFr* AsFr;
static mclBnFr* BsFr;
static mclBnFr* CsFr;

static mclBnFr* rsigma;
static mclBnFr* rsigmaInv;
static mclBnFr shift_fft;

static mpz_t* wM;

#include "../src/common/fourier.c"

int prover;
int cn;
int lro_constants_n;
int lro_const_total;
int wn;
int un;
int constant_n;

#include "../src/gro16/parser.c"

struct Trapdoor
{
    mpz_t alpha;
    mpz_t beta;
    mpz_t gamma;
    mpz_t delta;
    mpz_t x;
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

typedef struct
{
    mpz_t Ne;
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

struct mulExpResult
{
    mclBnG1 htdelta;
    mclBnG1 uwA1;
    mclBnG1 uwB1;
    mclBnG2 uwB2;
    mclBnG1 uwC1;
};

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
    proving_key pk;
    verifying_key vk;
} setup_keys;

#include "../src/common/utils.c"
#include "../src/gro16/prover.c"
#include "../src/gro16/qap.c"
#include "../src/gro16/setup.c"
#include "../src/gro16/verifier.c"

#endif