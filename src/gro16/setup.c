void setup(void *circuit, struct Trapdoor *t, struct Sigma1 *s1, struct Sigma2 *s2, mclBnGT *alphabetaT)
{
    mpz_t *A = (mpz_t*) malloc((M) * sizeof(mpz_t));
    mpz_t *B = (mpz_t*) malloc((M) * sizeof(mpz_t));
    mpz_t *C = (mpz_t*) malloc((M) * sizeof(mpz_t));

    mclBnG1 g;
    mclBnG2 h; // generators for G1 and G2
    mclBnG1_setStr(&g, GGEN, strlen(GGEN), 10);
    mclBnG2_setStr(&h, HGEN, strlen(HGEN), 10);

    mpz_init(t->alpha);
    mpz_init(t->beta);
    mpz_init(t->gamma);
    mpz_init(t->delta);
    mpz_init(t->x);

    mclBnFr rand;

    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&t->alpha, &rand);
    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&t->beta, &rand);
    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&t->gamma, &rand);
    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&t->delta, &rand);
    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&t->x, &rand);

    generateqap(circuit, A, B, C, *t);

    mpz_t factor, T;
    mpz_inits(factor, T, NULL);

    mpz_set_ui(factor, 1);
    mpz_powm(T, t->x, Ne, pPrime);
    mpz_sub(T, T, factor);

    // encrypt
    mclBnFr *frFactor = (mclBnFr*) malloc((M) * sizeof(mclBnFr));
    mpz_to_fr(&frFactor[0], &t->alpha);
    mclBnG1_mul(&s1->alpha, &g, &frFactor[0]);
    mpz_to_fr(&frFactor[0], &t->beta);
    mclBnG1_mul(&s1->beta, &g, &frFactor[0]);
    mpz_to_fr(&frFactor[0], &t->delta);
    mclBnG1_mul(&s1->delta, &g, &frFactor[0]);

    mpz_t invGamma, invDelta;
    mpz_inits(invGamma, invDelta, NULL);
    mpz_invert(invGamma, t->gamma, pPrime);
    mpz_invert(invDelta, t->delta, pPrime);

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_to_fr(&frFactor[i], &A[i]);
        mclBnG1_mul(&s1->A[i], &g, &frFactor[i]);
        mpz_to_fr(&frFactor[i], &B[i]);
        mclBnG1_mul(&s1->B[i], &g, &frFactor[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < nPublic; i++)
    {
        mpz_t f;
        mpz_init(f);
        // (t->beta * A[i] + t->alpha * B[i] + C[i]) * invGamma
        mpz_mul(f, t->beta, A[i]);
        mpz_mod(f, f, pPrime);
        mpz_addmul(f, t->alpha, B[i]);
        mpz_mod(f, f, pPrime);
        mpz_add(f, f, C[i]);

        mpz_mul(f, f, invGamma);
        mpz_mod(f, f, pPrime);
        mpz_to_fr(&frFactor[i], &f);
        mclBnG1_mul(&s1->vk[i], &g, &frFactor[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < M-nPublic; i++)
    {
        mpz_t f;
        mpz_init(f);
        // (t->beta * A[i] + t->alpha * B[i] + C[i]) * invDelta
        mpz_mul(f, t->beta, A[i+nPublic]);
        mpz_mod(f, f, pPrime);
        mpz_addmul(f, t->alpha, B[i+nPublic]);
        mpz_mod(f, f, pPrime);
        mpz_add(f, f, C[i+nPublic]);
        mpz_mul(f, f, invDelta);
        mpz_mod(f, f, pPrime);
        mpz_to_fr(&frFactor[i], &f);
        mclBnG1_mul(&s1->pk[i], &g, &frFactor[i]);
    }

    mpz_t zod, eT;
    mpz_inits(zod, eT, NULL);
    mpz_mul(zod, T, invDelta);
    mpz_mod(zod, zod, pPrime);

    mclBnFr *frFactor2 = (mclBnFr*) malloc((n) * sizeof(mclBnFr));

    mpz_to_fr(&frFactor2[0], &zod);
    mclBnG1_mul(&s1->xt[0], &g, &frFactor2[0]);
    mpz_set(eT, t->x);

    for (int i = 1; i < n; i++)
    {
        mpz_mul(factor, eT, zod);
        mpz_mod(factor, factor, pPrime);
        mpz_to_fr(&frFactor2[i], &factor);
        mclBnG1_mul(&s1->xt[i], &g, &frFactor2[i]);
        mpz_mul(eT, eT, t->x);
        mpz_mod(eT, eT, pPrime);
    }

    // encrypt
    mpz_to_fr(&frFactor[0], &t->beta);
    mclBnG2_mul(&s2->beta, &h, &frFactor[0]);
    mpz_to_fr(&frFactor[0], &t->gamma);
    mclBnG2_mul(&s2->gamma, &h, &frFactor[0]);
    mpz_to_fr(&frFactor[0], &t->delta);
    mclBnG2_mul(&s2->delta, &h, &frFactor[0]);

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_to_fr(&frFactor[i], &B[i]);
        mclBnG2_mul(&s2->B[i], &h, &frFactor[i]);
    }

    mclBn_pairing(alphabetaT, &s1->alpha, &s2->beta);
    
    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_clear(A[i]);
        mpz_clear(B[i]);
        mpz_clear(C[i]);
    }
}