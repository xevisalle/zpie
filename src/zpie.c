
void init_setup()
{
    M = 0;
    N = 0;
    nPublic = 0;

    mclBn_init(USEDCURVE, MCLBN_COMPILED_TIME_VAR);

    setParams = 1;
    init_circuit();
    setParams = 0;

    M+=100; // experimental

    uw = (mpz_t*) malloc((M) * sizeof(mpz_t));

    for (int i = 0; i < M; i++)
    {
        mpz_init2(uw[i], BITS);
    }
}

void perform_setup()
{
    struct stat st = {0};
    if (stat("data", &st) == -1) mkdir("data", 0700);

    struct Trapdoor t; // to be destroyed
    struct Sigma1 s1;
    struct Sigma2 s2;

    mpz_t kmul, factor, two;
    mpz_inits(kmul, Ne, factor, NULL);

    mpz_init(pPrime);
    mpz_set_str(pPrime, PRIMESTR, 10);
    mpz_sub_ui(factor, pPrime, 1);
    mpz_init_set_ui(two, 2);

    int i = 3;
    while ((mpz_cmp_ui(Ne, N) < 0) || (!mpz_divisible_p(factor, Ne)))
    {
        mpz_pow_ui(Ne, two, i);
        i++;
    }

    mpz_cdiv_q(kmul, factor, Ne);

    mpz_t base, w;
    mpz_init(w);
    mpz_init(base);
    mpz_init_set_ui(base, GROUPGEN); // multiplicative group generator  
    mpz_powm(w, base, kmul, pPrime);

    FILE *fpk, *fvk;
    fpk = fopen("data/provingkey.params", "w");
    fvk = fopen("data/verifyingkey.params", "w");

    mpz_out_str(fpk, 16, Ne);
    fprintf(fpk, "\n");

    n = mpz_get_ui(Ne);

    wM = (mpz_t*) malloc((n) * sizeof(mpz_t));
    s1.xt = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    s1.A = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.B = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.vk = (mclBnG1*) malloc((nPublic) * sizeof(mclBnG1));
    s1.pk = (mclBnG1*) malloc((M-nPublic) * sizeof(mclBnG1));
    s2.B = (mclBnG2*) malloc((M) * sizeof(mclBnG2));

    for (int i = 0; i < n; i++)
    {
        mpz_init(wM[i]);
        mpz_powm_ui(wM[i], w, i, pPrime);
        mpz_out_str(fpk, 16, wM[i]);
        fprintf(fpk, "\n");
    }

    mclBnGT alphabetaT;

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    setup(&t, &s1, &s2, &alphabetaT);

    fprintf(fpk, "%d\n", qapSize);

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(L[j][i]) fprintf(fpk, "%d\n%d\n%d\n", 1, j, i);
            if(R[j][i]) fprintf(fpk, "%d\n%d\n%d\n", 2, j, i);
            if(O[j][i]) fprintf(fpk, "%d\n%d\n%d\n", 3, j, i);
        }
    }

    char buff[2048];
    mclBnG1_getStr(buff, sizeof(buff), &s1.alpha, 16);
    fprintf(fpk, "%s\n", buff);
    mclBnG1_getStr(buff, sizeof(buff), &s1.beta, 16);
    fprintf(fpk, "%s\n", buff);
    mclBnG2_getStr(buff, sizeof(buff), &s2.beta, 16);
    fprintf(fpk, "%s\n", buff);
    mclBnG1_getStr(buff, sizeof(buff), &s1.delta, 16);
    fprintf(fpk, "%s\n", buff);
    mclBnG2_getStr(buff, sizeof(buff), &s2.delta, 16);
    fprintf(fpk, "%s\n", buff);

    for (int i = 0; i < M; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &s1.A[i], 16);
        fprintf(fpk, "%s\n", buff);
        mclBnG1_getStr(buff, sizeof(buff), &s1.B[i], 16);
        fprintf(fpk, "%s\n", buff);
        mclBnG2_getStr(buff, sizeof(buff), &s2.B[i], 16);
        fprintf(fpk, "%s\n", buff);
    }

    mclBnGT_getStr(buff, sizeof(buff), &alphabetaT, 10);
    fprintf(fvk, "%s\n", buff);

    mclBnG2_getStr(buff, sizeof(buff), &s2.gamma, 10);
    fprintf(fvk, "%s\n", buff);

    mclBnG2_getStr(buff, sizeof(buff), &s2.delta, 10);
    fprintf(fvk, "%s\n", buff);

    for (int i = 0; i < M-nPublic; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &s1.pk[i], 16);
        fprintf(fpk, "%s\n", buff);
    }

    for (int i = 0; i < nPublic; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &s1.vk[i], 10);
        fprintf(fvk, "%s\n", buff);
    }

    for (int i = 0; i < n; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &s1.xt[i], 16);
        fprintf(fpk, "%s\n", buff);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Setup generated successfully in", 1);
    if (bench) printf(" %fs\n", elapsed);
    
    fclose(fpk);
    fclose(fvk);
}

void init_prover()
{
    init_setup();

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    pk.A1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    pk.B1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    pk.B2 = (mclBnG2*) malloc((M) * sizeof(mclBnG2));
    pk.pk1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));

    if (bench) printf("  |--- Mode: Prove\n");

    char buff[2048];
    FILE *fpk;
    fpk = fopen("data/provingkey.params", "r");

    mpz_init(pPrime);
    mpz_set_str(pPrime, PRIMESTR, 10);
    fgets(buff, sizeof buff, fpk);
    mpz_init_set_str(Ne, buff, 16);

    n = mpz_get_ui(Ne);
    if (bench) printf("  |--- FFT constraints size : %d\n", n);

    wM = (mpz_t*) malloc((n) * sizeof(mpz_t));
    wMFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    AsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    BsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    CsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));

    for (int i = 0; i < n; i++)
    {
        fgets(buff, sizeof buff, fpk);
        mclBnFr_setStr(&wMFr[i], buff, strlen(buff), 16);
    }

    fgets(buff, sizeof buff, fpk);
    qapSize = atoi(buff);

    LRO = (int*) malloc((qapSize) * sizeof(int)); 

    for (int i = 0; i < qapSize; i++)
    {
        fgets(buff, sizeof buff, fpk);
        LRO[i] = atoi(buff);
    }

    fgets(buff, sizeof buff, fpk);
    mclBnG1_setStr(&pk.alpha1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG1_setStr(&pk.beta1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG2_setStr(&pk.beta2, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG1_setStr(&pk.delta1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG2_setStr(&pk.delta2, buff, strlen(buff), 16);

    for (int i = 0; i < M; i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&pk.A1[i], buff, strlen(buff), 16);
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&pk.B1[i], buff, strlen(buff), 16);
        fgets(buff,sizeof buff, fpk);
        mclBnG2_setStr(&pk.B2[i], buff, strlen(buff), 16);
    }

    for (int i = nPublic; i < M; i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&pk.pk1[i], buff, strlen(buff), 16);
    }

    pk.xt1 = (mclBnG1*) malloc((n) * sizeof(mclBnG1));

    for (int i = 0; i < n; i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&pk.xt1[i], buff, strlen(buff), 16);
    }

    fclose(fpk);   

    rsigma = (mpz_t*) malloc((n) * sizeof(mpz_t)); 
    rsigmaInv = (mpz_t*) malloc((n) * sizeof(mpz_t)); 

    mpz_t randNum;
    mpz_init(randNum);
    mpz_t factor;
    mpz_init_set_ui(factor, n);
    mpz_invert(factor, factor, pPrime);
    mpz_init(shift);

    mclBnFr rand;
    mclBnFr_setByCSPRNG(&rand);
    fr_to_mpz(&randNum, &rand);
    mpz_set(shift, randNum);

    mpz_init2(rsigma[0], BITS);
    mpz_init2(rsigmaInv[0], BITS);
    mpz_set_ui(rsigma[0], 1);
    mpz_invert(rsigmaInv[0], rsigma[0], pPrime);

    mclBnFr frFactor;
    mpz_to_fr(&frFactor, &rsigmaInv[0]);
    mclBnG1_mul(&pk.xt1[0], &pk.xt1[0], &frFactor);
    mpz_mul(rsigma[0], rsigma[0], factor);
    mpz_mod(rsigma[0], rsigma[0], pPrime);

    for (int i = 1; i < n; i++)
    {
        mclBnFr frFactorMulti;
        mpz_init2(rsigma[i], BITS);
        mpz_init2(rsigmaInv[i], BITS);
        mpz_powm_ui(rsigma[i], shift, i, pPrime);
        mpz_invert(rsigmaInv[i], rsigma[i], pPrime);

        mpz_to_fr(&frFactorMulti, &rsigmaInv[i]);
        mclBnG1_mul(&pk.xt1[i], &pk.xt1[i], &frFactorMulti);

        mpz_mul(rsigma[i], rsigma[i], factor);
        mpz_mod(rsigma[i], rsigma[i], pPrime);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("ZPiE started successfully in ", 1);
    if (bench) printf("%fs\n\n", elapsed);
}

proof generate_proof()
{
    uwn = 0;
    for (int i = 0; i < M; i++)
    {
        mpz_init(uw[i]);
    }
    
    proof p;

    if (bench) printf("--- Computing proof...\n");
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    prove(&p.piA, &p.piB2, &p.piC);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Proof generated successfully in ", 1);
    if (bench) printf("%fs\n", elapsed);

    p.uwProof = (mpz_t*) malloc((nPublic) * sizeof(mpz_t));

    for (int i = 0; i < nPublic; i++)
    {
        mpz_init(p.uwProof[i]);
        mpz_set(p.uwProof[i], uwProof[i]);
    }

    for (int i = 0; i < n; i++)
    {
        mclBnFr_clear(&AsFr[i]);
        mclBnFr_clear(&BsFr[i]);
        mclBnFr_clear(&CsFr[i]);
    }

    return p;
}

void store_proof(proof p)
{
    char buff[2048];
    FILE *fproof;
    fproof = fopen("data/proof.params", "w");

    for (int i = 0; i < nPublic; i++)
    {
        mpz_out_str(fproof, 10, p.uwProof[i]);
        fprintf(fproof, "\n");
    }

    mclBnG1_getStr(buff, sizeof(buff), &p.piA, 10);
    fprintf(fproof, "%s\n", buff);

    mclBnG2_getStr(buff, sizeof(buff), &p.piB2, 10);
    fprintf(fproof, "%s\n", buff);

    mclBnG1_getStr(buff, sizeof(buff), &p.piC, 10);
    fprintf(fproof, "%s\n", buff);

    fclose(fproof);
}

void init_verifier()
{
    init_setup();
    
    char buff[2048];
    FILE *fvk;
    fvk = fopen("data/verifyingkey.params", "r");

    fgets(buff, sizeof buff, fvk);
    mclBnGT_setStr(&vk.alphabetaT, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fvk);
    mclBnG2_setStr(&vk.gamma2, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fvk);
    mclBnG2_setStr(&vk.delta2, buff, strlen(buff), 10);

    vk.vk1 = (mclBnG1*) malloc((nPublic) * sizeof(mclBnG1));

    for (int i = 0; i < nPublic; i++)
    {
        fgets(buff, sizeof buff, fvk);
        mclBnG1_setStr(&vk.vk1[i], buff, strlen(buff), 10);
    }

    fclose(fvk);
}

proof read_proof()
{
    proof p;

    char buff[2048];
    FILE *fproof;
    fproof = fopen("data/proof.params", "r");

    p.uwProof = (mpz_t*) malloc((nPublic) * sizeof(mpz_t));

    for (int i = 0; i < nPublic; i++)
    {
        fgets(buff, sizeof buff, fproof);
        mpz_init(p.uwProof[i]);
        mpz_set_str(p.uwProof[i], buff, 10);
    }

    fgets(buff, sizeof buff, fproof);
    mclBnG1_setStr(&p.piA, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fproof);
    mclBnG2_setStr(&p.piB2, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fproof);
    mclBnG1_setStr(&p.piC, buff, strlen(buff), 10);

    fclose(fproof);

    return p;
}

int verify_proof(proof p)
{
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    int verified = verify(&p.piA, &p.piB2, &p.piC, p.uwProof);

    if (verified)
    {
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsed = (end.tv_sec - begin.tv_sec);
        elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
        log_success("Proof verified in ", 1);
        if (bench) printf("%fs\n", elapsed);
    }
    else
    {
        log_success("Proof cannot be verified\n", 0);
    }

    return verified;
}