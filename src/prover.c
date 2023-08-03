mpz_t *uwProof;

void h_coefficients()
{
    mclBnFr uwFr[M];

    #pragma omp parallel for
    for (int j = 0; j < M; j++)
    {
        mpz_to_fr(&uwFr[j], &uw[j]);
    }

    for (int j = 0; j < qapSize; j+=3)
    {
        switch (LRO[j])
        {
            case 1: mclBnFr_add(&AsFr[LRO[j+1]], &AsFr[LRO[j+1]], &uwFr[LRO[j+2]]); break;
            case 2: mclBnFr_add(&BsFr[LRO[j+1]], &BsFr[LRO[j+1]], &uwFr[LRO[j+2]]); break;
            case 3: mclBnFr_add(&CsFr[LRO[j+1]], &CsFr[LRO[j+1]], &uwFr[LRO[j+2]]); break;
        }
    }

    #pragma omp parallel num_threads(3)
    {
        switch (get_thread())
        {
            case 0: ifft_t(n, wMFr, AsFr); break;
            case 1: ifft_t(n, wMFr, BsFr); break;
            case 2: ifft_t(n, wMFr, CsFr); break;
            case 99:
                ifft_t(n, wMFr, AsFr);
                ifft_t(n, wMFr, BsFr);
                ifft_t(n, wMFr, CsFr);
                break;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {       
        mclBnFr_mul(&AsFr[i], &AsFr[i], &BsFr[i]);
        mclBnFr_sub(&AsFr[i], &AsFr[i], &CsFr[i]);
    }

    ifft(n, wMFr, AsFr);
}

void mul_exp(struct mulExpResult *result)
{
    int totTh = 16;
    uwProof = (mpz_t*) malloc((nPublic) * sizeof(mpz_t));

    for (int i = 0; i < nPublic; i++)
    {
        mpz_init(uwProof[i]);
        mpz_set(uwProof[i], uw[i]);
    }

    #ifdef MULTI
    mclBnG1 htdeltaTh[totTh];

    #pragma omp parallel num_threads(totTh)
    {
        mpz_t *exp[n/totTh];
        mclBnFr frFactor;
        int th = omp_get_thread_num();

        for (int i = th*(n/totTh); i < (th+1)*(n/totTh); i++)
        {
            fr_to_mpz(&wM[i], &AsFr[i]);
            exp[i-(th*(n/totTh))] = &wM[i];
        }

        bos_coster(exp, n/totTh, 1, &pk); 
        mpz_to_fr(&frFactor, exp[0]);
        mclBnG1_mul(&htdeltaTh[th], &pk.xt1[exp[0]-wM], &frFactor); 
    }

    mclBnG1_add(&result->htdelta, &htdeltaTh[0], &htdeltaTh[1]);

    for (int i = 2; i < totTh; i++)
    {
        mclBnG1_add(&result->htdelta, &result->htdelta, &htdeltaTh[i]);
    }

    mclBnG1 thA1[totTh];
    mclBnG1 thB1[totTh];
    mclBnG2 thB2[totTh];
    mclBnG1 thC1[totTh];

    #pragma omp parallel num_threads(totTh)
    {
        int sizeM;
        int end;
        int th = omp_get_thread_num();

        if (th == totTh-1)
        {
            sizeM = (M - (totTh*(M/totTh))) + (M/totTh);
            end = M;
        }
        else
        {
            sizeM = M/totTh;
            end = (th+1)*(M/totTh);
        }

        mclBnFr frFactor;
        mpz_t *exp[sizeM];

        int start = th*(M/totTh);

        for (int i = start; i < end; i++)
        {
            exp[i-start] = &uw[i];
        }

        bos_coster(exp, sizeM, 0, &pk);

        mpz_to_fr(&frFactor, exp[0]);
        mclBnG1_mul(&thA1[th], &pk.A1[exp[0]-uw], &frFactor);
        mclBnG1_mul(&thB1[th], &pk.B1[exp[0]-uw], &frFactor);
        mclBnG2_mul(&thB2[th], &pk.B2[exp[0]-uw], &frFactor);
        mclBnG1_mul(&thC1[th], &pk.pk1[exp[0]-uw], &frFactor);
    }

    mclBnG1_add(&result->uwA1, &thA1[0], &thA1[1]);
    mclBnG1_add(&result->uwB1, &thB1[0], &thB1[1]);
    mclBnG2_add(&result->uwB2, &thB2[0], &thB2[1]);
    mclBnG1_add(&result->uwC1, &thC1[0], &thC1[1]);

    for (int i = 2; i < totTh; i++)
    {
        mclBnG1_add(&result->uwA1, &result->uwA1, &thA1[i]);
        mclBnG1_add(&result->uwB1, &result->uwB1, &thB1[i]);
        mclBnG2_add(&result->uwB2, &result->uwB2, &thB2[i]);
        mclBnG1_add(&result->uwC1, &result->uwC1, &thC1[i]);
    }

    #else
    mpz_t *exp[n];
    mclBnFr frFactor;
    provingKey bpk;
    bpk.xt1 = (mclBnG1*) malloc((n) * sizeof(mclBnG1));

    for (int i = 0; i < n; i++)
    {
        fr_to_mpz(&wM[i], &AsFr[i]);
        exp[i] = &wM[i];
        mclBnG1_add(&bpk.xt1[i], &bpk.xt1[i], &pk.xt1[i]);
    }

    bos_coster(exp, n, 1, &bpk);
    mpz_to_fr(&frFactor, exp[0]);
    mclBnG1_mul(&result->htdelta, &bpk.xt1[exp[0]-wM], &frFactor);
    
    /*mpz_t *expM[M];
    bpk.A1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    bpk.B1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    bpk.B2 = (mclBnG2*) malloc((M) * sizeof(mclBnG2));
    bpk.pk1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        expM[i] = &uw[i];
        mclBnG1_add(&bpk.A1[i], &bpk.A1[i], &pk.A1[i]);
        mclBnG1_add(&bpk.B1[i], &bpk.B1[i], &pk.B1[i]);
        mclBnG2_add(&bpk.B2[i], &bpk.B2[i], &pk.B2[i]);
        mclBnG1_add(&bpk.pk1[i], &bpk.pk1[i], &pk.pk1[i]);
    }

    bos_coster(expM, M, 0, &bpk);

    mpz_to_fr(&frFactor, expM[0]);
    mclBnG1_mul(&result->uwA1, &bpk.A1[expM[0]-uw], &frFactor);
    mclBnG1_mul(&result->uwB1, &bpk.B1[expM[0]-uw], &frFactor);
    mclBnG2_mul(&result->uwB2, &bpk.B2[expM[0]-uw], &frFactor);
    mclBnG1_mul(&result->uwC1, &bpk.pk1[expM[0]-uw], &frFactor);*/

    // to be replaced ---->
    mclBnFr uwFactor[M];
    mclBnFr uwFactorPublic[M-nPublic];

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_to_fr(&uwFactor[i], &uw[i]);
        if(i >= nPublic) mpz_to_fr(&uwFactorPublic[i-nPublic], &uw[i]);
    }

    mclBnG1_mulVec(&result->uwA1, pk.A1, uwFactor, M);
    mclBnG1_mulVec(&result->uwB1, pk.B1, uwFactor, M);
    mclBnG2_mulVec(&result->uwB2, pk.B2, uwFactor, M);
    mclBnG1_mulVec(&result->uwC1, pk.pk1+nPublic, uwFactorPublic, M-nPublic);
    // <------ to be replaced
    #endif
}

void mcl_mul_exp(struct mulExpResult *result)
{
    mclBnFr uwFactor[M];
    mclBnFr uwFactorPublic[M-nPublic];

    uwProof = (mpz_t*) malloc((nPublic) * sizeof(mpz_t));

    for (int i = 0; i < nPublic; i++)
    {
        mpz_init(uwProof[i]);
        mpz_set(uwProof[i], uw[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_to_fr(&uwFactor[i], &uw[i]);
        if(i >= nPublic) mpz_to_fr(&uwFactorPublic[i-nPublic], &uw[i]);
    }

    #pragma omp parallel num_threads(5)
    {
        switch (get_thread())
        {
            case 0: mclBnG1_mulVec(&result->uwA1, pk.A1, uwFactor, M); break;
            case 1: mclBnG1_mulVec(&result->uwB1, pk.B1, uwFactor, M); break; 
            case 2: mclBnG2_mulVec(&result->uwB2, pk.B2, uwFactor, M); break;
            case 3: mclBnG1_mulVec(&result->uwC1, pk.pk1, uwFactorPublic, M-nPublic); break;
            case 4: mclBnG1_mulVec(&result->htdelta, pk.xt1, AsFr, n); break;
            case 99:
                mclBnG1_mulVec(&result->uwA1, pk.A1, uwFactor, M);
                mclBnG1_mulVec(&result->uwB1, pk.B1, uwFactor, M);
                mclBnG2_mulVec(&result->uwB2, pk.B2, uwFactor, M);
                mclBnG1_mulVec(&result->uwC1, pk.pk1+nPublic, uwFactorPublic, M-nPublic);
                mclBnG1_mulVec(&result->htdelta, pk.xt1, AsFr, n);
                break;
        }
    }
}

void naive_mul_exp(struct mulExpResult *result)
{
    mclBnFr frFactor[M];

    uwProof = (mpz_t*) malloc((nPublic) * sizeof(mpz_t));

    for (int i = 0; i < nPublic; i++)
    {
        mpz_init(uwProof[i]);
        mpz_set(uwProof[i], uw[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_to_fr(&frFactor[i], &uw[i]);
        // Auw = Auw + u[i] * s1.A[i];
        mclBnG1_mul(&pk.A1[i], &pk.A1[i], &frFactor[i]);
        // B1uw = B1uw + u[i] * s1.B[i];
        mclBnG1_mul(&pk.B1[i], &pk.B1[i], &frFactor[i]);
        // B2uw = B2uw + u[i] * s2.B[i];
        mclBnG2_mul(&pk.B2[i], &pk.B2[i], &frFactor[i]);
        // Cw = Cw + w[i] * s1.pk[i];
        if(i >= nPublic) mclBnG1_mul(&pk.pk1[i], &pk.pk1[i], &frFactor[i]);
    }

    mclBnG1_clear(&result->uwA1);
    mclBnG1_clear(&result->uwB1);
    mclBnG2_clear(&result->uwB2);
    mclBnG1_clear(&result->uwC1);

    for (int i = M; i--;)
    {
        mclBnG1_add(&result->uwA1, &result->uwA1, &pk.A1[i]);
        mclBnG1_add(&result->uwB1, &result->uwB1, &pk.B1[i]);
        mclBnG2_add(&result->uwB2, &result->uwB2, &pk.B2[i]);
        if(i >= nPublic) mclBnG1_add(&result->uwC1, &result->uwC1, &pk.pk1[i]);
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        mclBnG1_mul(&pk.xt1[i], &pk.xt1[i], &AsFr[i]);
    }

    mclBnG1_clear(&result->htdelta);

    for (int i = n; i--;)
    {
        mclBnG1_add(&result->htdelta, &result->htdelta, &pk.xt1[i]);
    }   
}

void prove(mclBnG1 *piA, mclBnG2 *piB2, mclBnG1 *piC)
{
    prover = 1;

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    init_circuit();

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    if (bench) printf("  |--- Circuit evaluation:  [%fs]\n", elapsed);

    clock_gettime(CLOCK_MONOTONIC, &begin);

    h_coefficients();

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    if (bench) printf("  |--- Compute h coefficients:  [%fs]\n", elapsed);

    clock_gettime(CLOCK_MONOTONIC, &begin);

    struct mulExpResult result;
    
    #ifdef AUTO_MULEXP
        if(M > 1000) mul_exp(&result);
        else mcl_mul_exp(&result);
    #elif BOSCOSTER_MULEXP
        mul_exp(&result);
    #elif NAIVE_MULEXP
        naive_mul_exp(&result);
    #elif MCL_MULEXP
        mcl_mul_exp(&result);
    #endif

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    if (bench) printf("  |--- G1, G2 multiexponentiations:  [%fs]\n", elapsed);
    #ifdef MULTI
    #elif BOSCOSTER_MULEXP
    if (bench) printf("     |--- Bos-Coster:  [%fs]\n", elapsedBosCoster);
    if (bench) printf("     |--- Heap sorting:  [%fs]\n", elapsedSort);
    #elif AUTO_MULEXP
    if (bench) printf("     |--- Bos-Coster:  [%fs]\n", elapsedBosCoster);
    if (bench) printf("     |--- Heap sorting:  [%fs]\n", elapsedSort);
    #endif

    log_message("Computing piA, piB1, piB2, piC, htdelta...");

    mclBnG1 piB1;
    mclBnFr r, s;

    mclBnFr_setByCSPRNG(&r);
    mclBnFr_setByCSPRNG(&s);

    // piA = s1.alpha + Auw + r * s1.delta;
    mclBnG1_mul(piA, &pk.delta1, &r);
    mclBnG1_add(piA, piA, &result.uwA1);
    mclBnG1_add(piA, piA, &pk.alpha1);
    // piB1 = s1.beta + B1uw + s * s1.delta;
    mclBnG1_mul(&piB1, &pk.delta1, &s);
    mclBnG1_add(&piB1, &piB1, &result.uwB1);
    mclBnG1_add(&piB1, &piB1, &pk.beta1);
    // piB2 = s2.beta + B2uw + s * s2.delta;
    mclBnG2_mul(piB2, &pk.delta2, &s);
    mclBnG2_add(piB2, piB2, &result.uwB2);
    mclBnG2_add(piB2, piB2, &pk.beta2);

    mclBnG1 factorG1;

    // piC = Cw + htdelta + piA*s + piB*r - r*s*s1.delta
    mclBnG1_mul(&factorG1, &pk.delta1, &r);
    mclBnG1_mul(&factorG1, &factorG1, &s);
    mclBnG1_mul(piC, &piB1, &r);
    mclBnG1_sub(&factorG1, piC, &factorG1);

    mclBnG1_mul(piC, piA, &s);
    mclBnG1_add(piC, piC, &factorG1);
    mclBnG1_add(piC, piC, &result.htdelta);
    mclBnG1_add(piC, piC, &result.uwC1);

    log_state(1);
    if (logs) 
    {
        char buff[2048];

        mclBnG1_getStr(buff, sizeof(buff), piA, 10);
        printf("piA = %s\n", buff);
        mclBnG1_getStr(buff, sizeof(buff), &piB1, 10);
        printf("piB1 = %s\n", buff);
        mclBnG2_getStr(buff, sizeof(buff), piB2, 10);
        printf("piB2 = %s\n", buff);
        mclBnG1_getStr(buff, sizeof(buff), piC, 10);
        printf("piC = %s\n", buff);
        mclBnG1_getStr(buff, sizeof(buff), &result.htdelta, 10);
        printf("htdelta = %s\n", buff);
    }
}