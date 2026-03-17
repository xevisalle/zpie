#include "../include/zpie.h"
#include <mcl/bn_c384_256.h>
#include <stdio.h>
#include <stdlib.h>

int logs;
int test_no_rand;

element one, oneNeg, c_mimc[91];
int bench;

int** L;
int** R;
int** O;

mclBnFr* LRO_constants;

static mclBnFr* AsFr;
static mclBnFr* BsFr;
static mclBnFr* CsFr;

static mclBnFr* rsigma;
static mclBnFr* rsigmaInv;
static mclBnFr shift_fft;

static mclBnFr* wM;

int logN;
int Nb;
int Mc;

int M;
int N;
int nPublic;
int nConst;
int setParams;

mclBnFr* uw;

#include "../src/common/fourier.c"

int prover;
int cn;
int lro_constants_n;
int lro_const_total;
int wn;
int un;
int constant_n;

#include "../src/common/utils.c"
#include "../src/gro16/parser.c"
#include "../src/gro16/prover.c"
#include "../src/gro16/qap.c"
#include "../src/gro16/setup.c"
#include "../src/gro16/verifier.c"

setup_keys perform_setup(void* circuit)
{
    init_setup(circuit);

    struct Trapdoor t; // to be destroyed

    struct Sigma1 s1;
    struct Sigma2 s2;
    mclBnGT alphabetaT;

    // Find smallest power of 2 >= N that divides (p - 1)
    // Since p - 1 = r * 2^s for BN/BLS curves, any power of 2 up to 2^s works.
    // We compute this with plain integers since Ne is always small.
    mclBnFr frPrime;
    mclBnFr_setStr(&frPrime, PRIMESTR, strlen(PRIMESTR), 10);

    int Ne = 8; // start at 2^3
    while (Ne < N)
    {
        Ne <<= 1;
    }
    // Ne is now a power of 2 >= N.
    // For BN128/BLS12_381, (p-1) is divisible by large powers of 2, so this is fine.

    // Compute root of unity: w = g^((p-1)/Ne) mod p
    // Using mclBnFr: factor = (p-1)/Ne, w = GROUPGEN^factor
    mclBnFr frOne, frNe, factor, w;
    mclBnFr_setInt(&frOne, 1);
    mclBnFr_setInt(&frNe, Ne);

    // (p - 1) in Fr: since Fr is mod p, p mod p = 0, so p-1 = -1 in Fr
    mclBnFr pMinus1;
    mclBnFr_neg(&pMinus1, &frOne); // -1 mod p = p-1

    mclBnFr_div(&factor, &pMinus1, &frNe); // (p-1)/Ne

    mclBnFr base;
    mclBnFr_setInt(&base, GROUPGEN);
    mclBnFr_pow(&w, &base, &factor);

    int n = Ne;

    setup_keys keys;
    keys.pk.Ne = Ne;

    keys.pk.LRO_constants = (mclBnFr*) malloc((lro_const_total) * sizeof(mclBnFr));
    keys.pk.wM = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

    wM = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    s1.xt = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    s1.A = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.B = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.vk = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));
    s1.pk = (mclBnG1*) malloc((M - (nPublic + nConst)) * sizeof(mclBnG1));
    s2.B = (mclBnG2*) malloc((M) * sizeof(mclBnG2));

    // wM[i] = w^i
    mclBnFr_setInt(&wM[0], 1);
    for (int i = 1; i < n; i++)
    {
        mclBnFr_mul(&wM[i], &wM[i - 1], &w);
    }
    for (int i = 0; i < n; i++)
    {
        keys.pk.wM[i] = wM[i];
    }

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    setup(circuit, &t, &s1, &s2, &alphabetaT, &keys.pk.qap_size, keys.pk.Ne);

    keys.pk.LRO = (int*) malloc((keys.pk.qap_size) * sizeof(int));

    for (int i = 0; i < lro_const_total; i++)
    {
        keys.pk.LRO_constants[i] = LRO_constants[i];
    }

    int it = 0;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (L[j][i] != 0)
            {
                keys.pk.LRO[it + 1] = j;
                keys.pk.LRO[it + 2] = i;

                if (L[j][i] != 1)
                {
                    keys.pk.LRO[it] = 10;
                    keys.pk.LRO[it + 3] = L[j][i];
                    it += 4;
                }
                else
                {
                    keys.pk.LRO[it] = 1;
                    it += 3;
                }
            }
            if (R[j][i] != 0)
            {
                keys.pk.LRO[it + 1] = j;
                keys.pk.LRO[it + 2] = i;

                if (R[j][i] != 1)
                {
                    keys.pk.LRO[it] = 20;
                    keys.pk.LRO[it + 3] = R[j][i];
                    it += 4;
                }
                else
                {
                    keys.pk.LRO[it] = 2;
                    it += 3;
                }
            }
            if (O[j][i])
            {
                keys.pk.LRO[it] = 3;
                keys.pk.LRO[it + 1] = j;
                keys.pk.LRO[it + 2] = i;
                it += 3;
            }
        }
    }

    keys.pk.alpha1 = s1.alpha;
    keys.pk.beta1 = s1.beta;
    keys.pk.beta2 = s2.beta;
    keys.pk.delta1 = s1.delta;
    keys.pk.delta2 = s2.delta;
    keys.pk.A1 = s1.A;
    keys.pk.B1 = s1.B;
    keys.pk.B2 = s2.B;
    keys.pk.pk1 = s1.pk;
    keys.pk.xt1 = s1.xt;
    keys.pk.xt1_rand = (mclBnG1*) malloc((n) * sizeof(mclBnG1));

    keys.vk.alphabetaT = alphabetaT;
    keys.vk.gamma2 = s2.gamma;
    keys.vk.delta2 = s2.delta;

    keys.vk.constants = (mclBnFr*) malloc((nConst) * sizeof(mclBnFr));

    for (int i = 0; i < (nConst); i++)
    {
        mclBnFr_neg(&keys.vk.constants[i], &uw[i]);
        mclBnFr_neg(&keys.vk.constants[i], &keys.vk.constants[i]);
    }

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        keys.vk.vk1[i] = s1.vk[i];
    }

    free(s1.vk);

    for (int i = 0; i < N; i++)
    {
        free(L[i]);
        free(R[i]);
        free(O[i]);
    }
    free(L);
    free(R);
    free(O);
    L = NULL;
    R = NULL;
    O = NULL;

    free(wM);
    wM = NULL;

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Setup generated successfully in", 1);
    if (bench)
        printf(" %fs\n", elapsed);

    return keys;
}

void serialize_pk(proving_key* pk)
{
    FILE* fpk;
    fpk = fopen("data/provingkey.params", "w");

    int n = pk->Ne;

    int buff_pk_size = SIZE_FR * (n + lro_const_total) + SIZE_G2 * (2 + M) +
                       SIZE_G1 * (M - (nPublic + nConst) + 3 + n + 2 * M);
    char buff_pk[buff_pk_size];

    fwrite(&pk->Ne, sizeof(int), 1, fpk);
    fwrite(&pk->qap_size, sizeof(int), 1, fpk);

    for (int i = 0; i < pk->qap_size; i++)
    {
        fwrite(&pk->LRO[i], sizeof(int), 1, fpk);
    }

    int size = 0;

    for (int i = 0; i < lro_const_total; i++)
    {
        size += mclBnFr_serialize(buff_pk + size, SIZE_FR, &pk->LRO_constants[i]);
    }

    for (int i = 0; i < n; i++)
    {
        size += mclBnFr_serialize(buff_pk + size, SIZE_FR, &pk->wM[i]);
    }

    size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->alpha1);
    size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->beta1);
    size += mclBnG2_serialize(buff_pk + size, SIZE_G2, &pk->beta2);
    size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->delta1);
    size += mclBnG2_serialize(buff_pk + size, SIZE_G2, &pk->delta2);

    for (int i = 0; i < M; i++)
    {
        size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->A1[i]);
        size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->B1[i]);
        size += mclBnG2_serialize(buff_pk + size, SIZE_G2, &pk->B2[i]);
    }

    for (int i = 0; i < M - (nPublic + nConst); i++)
    {
        size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->pk1[i]);
    }

    for (int i = 0; i < n; i++)
    {
        size += mclBnG1_serialize(buff_pk + size, SIZE_G1, &pk->xt1[i]);
    }

    fwrite(buff_pk, 1, size, fpk);
    fclose(fpk);
}

void serialize_vk(verifying_key* vk)
{
    FILE* fvk;
    fvk = fopen("data/verifyingkey.params", "w");

    int buff_vk_size = SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst) + SIZE_FR * nConst;
    char buff_vk[buff_vk_size];

    int size = 0;

    for (int i = 0; i < nConst; i++)
    {
        size += mclBnFr_serialize(buff_vk + size, SIZE_FR, &vk->constants[i]);
    }

    size += mclBnGT_serialize(buff_vk + size, SIZE_GT, &vk->alphabetaT);
    size += mclBnG2_serialize(buff_vk + size, SIZE_G2, &vk->gamma2);
    size += mclBnG2_serialize(buff_vk + size, SIZE_G2, &vk->delta2);

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        size += mclBnG1_serialize(buff_vk + size, SIZE_G1, &vk->vk1[i]);
    }

    fwrite(buff_vk, 1, size, fvk);
    fclose(fvk);
}

void store_setup(setup_keys* keys)
{
    struct stat st = {0};
    if (stat("data", &st) == -1)
        mkdir("data", 0700);

    serialize_pk(&keys->pk);
    serialize_vk(&keys->vk);
}

setup_keys read_setup(void* circuit)
{
    init_setup(circuit);

    FILE *fpk, *fvk;

    fpk = fopen("data/provingkey.params", "r");
    fvk = fopen("data/verifyingkey.params", "r");

    setup_keys keys;

    size_t __attribute__((unused)) _nr;
    _nr = fread(&keys.pk.Ne, sizeof(int), 1, fpk);

    int n = keys.pk.Ne;
    int buff_pk_size = SIZE_FR * (n + lro_const_total) + SIZE_G2 * (2 + M) +
                       SIZE_G1 * (M - (nPublic + nConst) + 3 + n + 2 * M);
    char buff_pk[buff_pk_size];

    keys.pk.wM = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    keys.vk.constants = (mclBnFr*) malloc(((nConst)) * sizeof(mclBnFr));

    keys.pk.xt1 = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.xt1_rand = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.A1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.B1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.pk1 = (mclBnG1*) malloc((M - (nPublic + nConst)) * sizeof(mclBnG1));
    keys.pk.B2 = (mclBnG2*) malloc((M) * sizeof(mclBnG2));
    keys.pk.LRO_constants = (mclBnFr*) malloc((lro_const_total) * sizeof(mclBnFr));

    _nr = fread(&keys.pk.qap_size, sizeof(int), 1, fpk);
    keys.pk.LRO = (int*) malloc((keys.pk.qap_size) * sizeof(int));

    for (int i = 0; i < keys.pk.qap_size; i++)
    {
        _nr = fread(&keys.pk.LRO[i], sizeof(int), 1, fpk);
    }

    int size = 0;
    _nr = fread(buff_pk, 1, buff_pk_size, fpk);

    for (int i = 0; i < lro_const_total; i++)
    {
        size += mclBnFr_deserialize(&keys.pk.LRO_constants[i], buff_pk + size, SIZE_FR);
    }

    for (int i = 0; i < n; i++)
    {
        size += mclBnFr_deserialize(&keys.pk.wM[i], buff_pk + size, SIZE_FR);
    }

    size += mclBnG1_deserialize(&keys.pk.alpha1, buff_pk + size, SIZE_G1);
    size += mclBnG1_deserialize(&keys.pk.beta1, buff_pk + size, SIZE_G1);
    size += mclBnG2_deserialize(&keys.pk.beta2, buff_pk + size, SIZE_G2);
    size += mclBnG1_deserialize(&keys.pk.delta1, buff_pk + size, SIZE_G1);
    size += mclBnG2_deserialize(&keys.pk.delta2, buff_pk + size, SIZE_G2);

    for (int i = 0; i < M; i++)
    {
        size += mclBnG1_deserialize(&keys.pk.A1[i], buff_pk + size, SIZE_G1);
        size += mclBnG1_deserialize(&keys.pk.B1[i], buff_pk + size, SIZE_G1);
        size += mclBnG2_deserialize(&keys.pk.B2[i], buff_pk + size, SIZE_G2);
    }

    for (int i = 0; i < M - (nPublic + nConst); i++)
    {
        size += mclBnG1_deserialize(&keys.pk.pk1[i], buff_pk + size, SIZE_G1);
    }

    for (int i = 0; i < n; i++)
    {
        size += mclBnG1_deserialize(&keys.pk.xt1[i], buff_pk + size, SIZE_G1);
    }

    char buff_vk[SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst) + SIZE_FR * nConst];
    size = 0;

    _nr = fread(buff_vk, 1, SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst) + SIZE_FR * nConst, fvk);

    for (int i = 0; i < nConst; i++)
    {
        size += mclBnFr_deserialize(&keys.vk.constants[i], buff_vk + size, SIZE_FR);
    }

    size += mclBnGT_deserialize(&keys.vk.alphabetaT, buff_vk + size, SIZE_GT);
    size += mclBnG2_deserialize(&keys.vk.gamma2, buff_vk + size, SIZE_G2);
    size += mclBnG2_deserialize(&keys.vk.delta2, buff_vk + size, SIZE_G2);

    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        size += mclBnG1_deserialize(&keys.vk.vk1[i], buff_vk + size, SIZE_G1);
    }

    fclose(fpk);
    fclose(fvk);

    return keys;
}

proof generate_proof(void* circuit, proving_key* pk)
{
    init_prover(circuit, pk);

    wn = nPublic + nConst;
    un = nConst;
    constant_n = 0;
    for (int i = 0; i < M; i++)
    {
        mclBnFr_clear(&uw[i]);
    }

    int n = pk->Ne;

    proof p;

    p.uwProof = (mclBnFr*) malloc((nPublic) * sizeof(mclBnFr));

    if (bench)
        printf("--- Computing proof...\n");
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    prove(circuit, &p.piA, &p.piB2, &p.piC, p.uwProof, pk);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Proof generated successfully in ", 1);
    if (bench)
        printf("%fs\n", elapsed);

    free(AsFr);
    free(BsFr);
    free(CsFr);
    free(rsigma);
    free(rsigmaInv);
    AsFr = NULL;
    BsFr = NULL;
    CsFr = NULL;
    rsigma = NULL;
    rsigmaInv = NULL;

    return p;
}

void store_proof(proof* p)
{
    char buff[2048];
    FILE* fproof;
    fproof = fopen("data/proof.params", "w");

    int size = 0;

    for (int i = 0; i < (nPublic); i++)
    {
        size += mclBnFr_serialize(buff + size, SIZE_FR, &p->uwProof[i]);
    }

    size += mclBnG1_serialize(buff + size, SIZE_G1, &p->piA);
    size += mclBnG2_serialize(buff + size, SIZE_G2, &p->piB2);
    size += mclBnG1_serialize(buff + size, SIZE_G1, &p->piC);

    fwrite(buff, 1, size, fproof);
    fclose(fproof);
}

proof read_proof()
{
    proof p;

    char buff[2048];
    FILE* fproof;
    fproof = fopen("data/proof.params", "r");

    p.uwProof = (mclBnFr*) malloc((nPublic) * sizeof(mclBnFr));

    int size = 0;

    size_t __attribute__((unused)) _nr;
    _nr = fread(buff, 1, (SIZE_FR * nPublic) + SIZE_G1 + SIZE_G2 + SIZE_G1, fproof);

    for (int i = 0; i < (nPublic); i++)
    {
        size += mclBnFr_deserialize(&p.uwProof[i], buff + size, SIZE_FR);
    }

    size += mclBnG1_deserialize(&p.piA, buff + size, SIZE_G1);
    size += mclBnG2_deserialize(&p.piB2, buff + size, SIZE_G2);
    size += mclBnG1_deserialize(&p.piC, buff + size, SIZE_G1);

    fclose(fproof);

    return p;
}

int verify_proof(void* circuit, proof* p, verifying_key* vk)
{
    init_setup(circuit);

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    int verified = verify(p, vk);

    if (verified)
    {
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsed = (end.tv_sec - begin.tv_sec);
        elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
        log_success("Proof verified in ", 1);
        if (bench)
            printf("%fs\n", elapsed);
    }
    else
    {
        log_success("Proof cannot be verified\n", 0);
    }

    return verified;
}
