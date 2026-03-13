void setup(void* circuit, struct Trapdoor* t, struct Sigma1* s1, struct Sigma2* s2,
           mclBnGT* alphabetaT, int* qap_size, int Ne)
{
    mclBnFr* A = (mclBnFr*) malloc((M) * sizeof(mclBnFr));
    mclBnFr* B = (mclBnFr*) malloc((M) * sizeof(mclBnFr));
    mclBnFr* C = (mclBnFr*) malloc((M) * sizeof(mclBnFr));

    mclBnG1 g;
    mclBnG2 h; // generators for G1 and G2
    mclBnG1_setStr(&g, GGEN, strlen(GGEN), 10);
    mclBnG2_setStr(&h, HGEN, strlen(HGEN), 10);

    generate_random_scalar(&t->alpha);
    generate_random_scalar(&t->beta);
    generate_random_scalar(&t->gamma);
    generate_random_scalar(&t->delta);
    generate_random_scalar(&t->x);

    generateqap(circuit, A, B, C, *t, qap_size, Ne);

    mclBnFr frNe, T, frOne;
    mclBnFr_setInt(&frNe, Ne);
    mclBnFr_setInt(&frOne, 1);
    mclBnFr_pow(&T, &t->x, &frNe);
    mclBnFr_sub(&T, &T, &frOne);

    // encrypt
    mclBnG1_mul(&s1->alpha, &g, &t->alpha);
    mclBnG1_mul(&s1->beta, &g, &t->beta);
    mclBnG1_mul(&s1->delta, &g, &t->delta);

    mclBnFr invGamma, invDelta;
    mclBnFr_inv(&invGamma, &t->gamma);
    mclBnFr_inv(&invDelta, &t->delta);

#pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mclBnG1_mul(&s1->A[i], &g, &A[i]);
        mclBnG1_mul(&s1->B[i], &g, &B[i]);
    }

#pragma omp parallel for
    for (int i = 0; i < (nPublic + nConst); i++)
    {
        mclBnFr f, tmp;
        // (t->beta * A[i] + t->alpha * B[i] + C[i]) * invGamma
        mclBnFr_mul(&f, &t->beta, &A[i]);
        mclBnFr_mul(&tmp, &t->alpha, &B[i]);
        mclBnFr_add(&f, &f, &tmp);
        mclBnFr_add(&f, &f, &C[i]);
        mclBnFr_mul(&f, &f, &invGamma);
        mclBnG1_mul(&s1->vk[i], &g, &f);
    }

#pragma omp parallel for
    for (int i = 0; i < M - (nPublic + nConst); i++)
    {
        mclBnFr f, tmp;
        // (t->beta * A[i] + t->alpha * B[i] + C[i]) * invDelta
        mclBnFr_mul(&f, &t->beta, &A[i + (nPublic + nConst)]);
        mclBnFr_mul(&tmp, &t->alpha, &B[i + (nPublic + nConst)]);
        mclBnFr_add(&f, &f, &tmp);
        mclBnFr_add(&f, &f, &C[i + (nPublic + nConst)]);
        mclBnFr_mul(&f, &f, &invDelta);
        mclBnG1_mul(&s1->pk[i], &g, &f);
    }

    mclBnFr zod, eT;
    mclBnFr_mul(&zod, &T, &invDelta);

    int n = Ne;

    mclBnG1_mul(&s1->xt[0], &g, &zod);
    eT = t->x;

    for (int i = 1; i < n; i++)
    {
        mclBnFr factor;
        mclBnFr_mul(&factor, &eT, &zod);
        mclBnG1_mul(&s1->xt[i], &g, &factor);
        mclBnFr_mul(&eT, &eT, &t->x);
    }

    // encrypt
    mclBnG2_mul(&s2->beta, &h, &t->beta);
    mclBnG2_mul(&s2->gamma, &h, &t->gamma);
    mclBnG2_mul(&s2->delta, &h, &t->delta);

#pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mclBnG2_mul(&s2->B[i], &h, &B[i]);
    }

    mclBn_pairing(alphabetaT, &s1->alpha, &s2->beta);

    free(A);
    free(B);
    free(C);
}