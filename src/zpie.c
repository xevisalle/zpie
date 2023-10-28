
setup_keys perform_setup(void *circuit)
{
    init_setup(circuit);
    
    struct Trapdoor t; // to be destroyed

    struct Sigma1 s1;
    struct Sigma2 s2;
    mclBnGT alphabetaT;

    mpz_t kmul, factor, two, Ne;
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

    int n = mpz_get_ui(Ne);

    setup_keys provk;
    mpz_init_set(provk.pk.Ne, Ne);

    provk.pk.wMFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    provk.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

    wM = (mpz_t*) malloc((n) * sizeof(mpz_t));
    s1.xt = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    s1.A = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.B = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    s1.vk = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));
    s1.pk = (mclBnG1*) malloc((M-(nPublic + nConst)) * sizeof(mclBnG1));
    s2.B = (mclBnG2*) malloc((M) * sizeof(mclBnG2));

    for (int i = 0; i < n; i++)
    {
        mpz_init(wM[i]);
        mpz_powm_ui(wM[i], w, i, pPrime);
        mpz_to_fr(&provk.pk.wMFr[i], &wM[i]);
    }

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    setup(circuit, &t, &s1, &s2, &alphabetaT, &provk.pk.qap_size, &provk.pk.Ne);

    provk.pk.LRO = (int*) malloc((provk.pk.qap_size) * sizeof(int));

    int it = 0;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(L[j][i] != 0) 
            {
                provk.pk.LRO[it+1] = j;
                provk.pk.LRO[it+2] = i;

                if (L[j][i] != 1)
                {
                    provk.pk.LRO[it] = 10;
                    provk.pk.LRO[it+3] = L[j][i];
                    it+=4;
                }
                else
                {
                    provk.pk.LRO[it] = 1;
                    it+=3;
                }
            }
            if(R[j][i] != 0) 
            {
                provk.pk.LRO[it+1] = j;
                provk.pk.LRO[it+2] = i;

                if (R[j][i] != 1)
                {
                    provk.pk.LRO[it] = 20;
                    provk.pk.LRO[it+3] = R[j][i];
                    it+=4;
                }
                else 
                {
                    provk.pk.LRO[it] = 2;
                    it+=3;
                }
            }
            if(O[j][i]) 
            {
                provk.pk.LRO[it] = 3;
                provk.pk.LRO[it+1] = j;
                provk.pk.LRO[it+2] = i;
                it+=3;
            }
        }
    }

    provk.pk.alpha1 = s1.alpha;
    provk.pk.beta1 = s1.beta;
    provk.pk.beta2 = s2.beta;
    provk.pk.delta1 = s1.delta;
    provk.pk.delta2 = s2.delta;
    provk.pk.A1 = s1.A;
    provk.pk.B1 = s1.B;
    provk.pk.B2 = s2.B;
    provk.pk.pk1 = s1.pk;
    provk.pk.xt1 = s1.xt;
    provk.pk.xt1_rand = (mclBnG1*) malloc((n) * sizeof(mclBnG1));

    provk.vk.alphabetaT = alphabetaT;
    provk.vk.gamma2 = s2.gamma;
    provk.vk.delta2 = s2.delta;


    for (int i = 0; i < (nPublic + nConst); i++)
    {
        provk.vk.vk1[i] = s1.vk[i];
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Setup generated successfully in", 1);
    if (bench) printf(" %fs\n", elapsed);

    return provk;
}

char* serialize_pk(proving_key *pk)
{
    int n = mpz_get_ui(pk->Ne);

    char *pk_bytes;
    pk_bytes = (char *) malloc(1024 * n * sizeof(char));

    for (int i = 0; i < 1024 * n * sizeof(char); i++)
    {
        pk_bytes[i] = 0;
    }

    char buff[2048];
    mpz_get_str(buff, 16, pk->Ne);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");

    for (int i = 0; i < n; i++)
    {
        mclBnFr_getStr(buff, sizeof(buff), &pk->wMFr[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
    }

    sprintf(buff, "%d", pk->qap_size);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");

    for (int i = 0; i < pk->qap_size; i++)
    {
        sprintf(buff, "%d", pk->LRO[i]);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
    }

    mclBnG1_getStr(buff, sizeof(buff), &pk->alpha1, 16);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");
    mclBnG1_getStr(buff, sizeof(buff), &pk->beta1, 16);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");
    mclBnG2_getStr(buff, sizeof(buff), &pk->beta2, 16);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");
    mclBnG1_getStr(buff, sizeof(buff), &pk->delta1, 16);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");
    mclBnG2_getStr(buff, sizeof(buff), &pk->delta2, 16);
    strcat(pk_bytes, buff);
    strcat(pk_bytes, "\n");

    for (int i = 0; i < M; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &pk->A1[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
        mclBnG1_getStr(buff, sizeof(buff), &pk->B1[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
        mclBnG2_getStr(buff, sizeof(buff), &pk->B2[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
    }

    for (int i = 0; i < M-(nPublic + nConst); i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &pk->pk1[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
    }

    for (int i = 0; i < n; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &pk->xt1[i], 16);
        strcat(pk_bytes, buff);
        strcat(pk_bytes, "\n");
    }

    return pk_bytes;
}

char* serialize_vk(verifying_key *vk)
{
    char *vk_bytes;
    vk_bytes = (char *) malloc(1024 * (nPublic + nConst) * sizeof(char));

    for (int i = 0; i < 1024 * (nPublic + nConst) * sizeof(char); i++)
    {
        vk_bytes[i] = 0;
    }

    char buff[2048];
    mclBnGT_getStr(buff, sizeof(buff), &vk->alphabetaT, 10);
    strcat(vk_bytes, buff);
    strcat(vk_bytes, "\n");

    mclBnG2_getStr(buff, sizeof(buff), &vk->gamma2, 10);
    strcat(vk_bytes, buff);
    strcat(vk_bytes, "\n");

    mclBnG2_getStr(buff, sizeof(buff), &vk->delta2, 10);
    strcat(vk_bytes, buff);
    strcat(vk_bytes, "\n");

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &vk->vk1[i], 10);
        strcat(vk_bytes, buff);
        strcat(vk_bytes, "\n");
    }

    return vk_bytes;
}

void store_setup(setup_keys keys)
{
    struct stat st = {0};
    if (stat("data", &st) == -1) mkdir("data", 0700);

    FILE *fpk, *fvk;
    fpk = fopen("data/provingkey.params", "w");
    fvk = fopen("data/verifyingkey.params", "w");

    char *pk_bytes = serialize_pk(&keys.pk);
    char *vk_bytes = serialize_vk(&keys.vk);

    fprintf(fpk, "%s", pk_bytes);
    fprintf(fvk, "%s", vk_bytes);
    
    fclose(fpk);
    fclose(fvk);
}

setup_keys read_setup(void *circuit)
{
    init_setup(circuit);

    char buff[2048];
    FILE *fpk, *fvk;
    
    fpk = fopen("data/provingkey.params", "r");
    fvk = fopen("data/verifyingkey.params", "r");

    setup_keys keys;

    fgets(buff, sizeof buff, fpk);
    mpz_init_set_str(keys.pk.Ne, buff, 16);

    int n = mpz_get_ui(keys.pk.Ne);
    
    keys.pk.wMFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

    keys.pk.xt1 = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.xt1_rand = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.A1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.B1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.pk1 = (mclBnG1*) malloc((M-(nPublic + nConst)) * sizeof(mclBnG1));
    keys.pk.B2 = (mclBnG2*) malloc((M) * sizeof(mclBnG2));

    for (int i = 0; i < n; i++)
    {
        fgets(buff, sizeof buff, fpk);
        mclBnFr_setStr(&keys.pk.wMFr[i], buff, strlen(buff), 16);
    }

    fgets(buff, sizeof buff, fpk);
    keys.pk.qap_size = atoi(buff);

    keys.pk.LRO = (int*) malloc((keys.pk.qap_size) * sizeof(int)); 

    for (int i = 0; i < keys.pk.qap_size; i++)
    {
        fgets(buff, sizeof buff, fpk);
        keys.pk.LRO[i] = atoi(buff);
    }

    fgets(buff, sizeof buff, fpk);
    mclBnG1_setStr(&keys.pk.alpha1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG1_setStr(&keys.pk.beta1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG2_setStr(&keys.pk.beta2, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG1_setStr(&keys.pk.delta1, buff, strlen(buff), 16);
    fgets(buff,sizeof buff, fpk);
    mclBnG2_setStr(&keys.pk.delta2, buff, strlen(buff), 16);

    for (int i = 0; i < M; i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&keys.pk.A1[i], buff, strlen(buff), 16);
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&keys.pk.B1[i], buff, strlen(buff), 16);
        fgets(buff,sizeof buff, fpk);
        mclBnG2_setStr(&keys.pk.B2[i], buff, strlen(buff), 16);
    }

    for (int i = 0; i < M-(nPublic + nConst); i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&keys.pk.pk1[i], buff, strlen(buff), 16);
    }

    for (int i = 0; i < n; i++)
    {
        fgets(buff,sizeof buff, fpk);
        mclBnG1_setStr(&keys.pk.xt1[i], buff, strlen(buff), 16);
    }
    
    fgets(buff, sizeof buff, fvk);
    mclBnGT_setStr(&keys.vk.alphabetaT, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fvk);
    mclBnG2_setStr(&keys.vk.gamma2, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fvk);
    mclBnG2_setStr(&keys.vk.delta2, buff, strlen(buff), 10);

    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        fgets(buff, sizeof buff, fvk);
        mclBnG1_setStr(&keys.vk.vk1[i], buff, strlen(buff), 10);
    }

    fclose(fpk); 
    fclose(fvk);

    return keys;
}

proof generate_proof(void *circuit, proving_key pk)
{
    init_prover(circuit, pk);

    wn = nPublic + nConst;
    un = nConst;
    constant_n = 0;
    for (int i = 0; i < M; i++)
    {
        mpz_init(uw[i]);
    }

    int n = mpz_get_ui(pk.Ne);
    wM = (mpz_t*) malloc((n) * sizeof(mpz_t));
    
    proof p;

    p.uwProof = (mpz_t*) malloc((nPublic+nConst) * sizeof(mpz_t));

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        mpz_init(p.uwProof[i]);
    }

    if (bench) printf("--- Computing proof...\n");
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    prove(circuit, &p.piA, &p.piB2, &p.piC, p.uwProof, pk);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Proof generated successfully in ", 1);
    if (bench) printf("%fs\n", elapsed);

    for (int i = 0; i < n; i++)
    {
        mclBnFr_clear(&AsFr[i]);
        mclBnFr_clear(&BsFr[i]);
        mclBnFr_clear(&CsFr[i]);

        mpz_clear(rsigma[i]);
        mpz_clear(rsigmaInv[i]);
    }

    for (int i = 0; i < M; i++)
    {
        mpz_clear(uw[i]);
    }

    mpz_clear(shift);

    return p;
}

void store_proof(proof p)
{
    char buff[2048];
    FILE *fproof;
    fproof = fopen("data/proof.params", "w");

    for (int i = 0; i < (nPublic + nConst); i++)
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

proof read_proof()
{
    proof p;

    char buff[2048];
    FILE *fproof;
    fproof = fopen("data/proof.params", "r");

    p.uwProof = (mpz_t*) malloc((nPublic+nConst) * sizeof(mpz_t));

    for (int i = 0; i < (nPublic + nConst); i++)
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

int verify_proof(void *circuit, proof p, verifying_key vk)
{
    init_setup(circuit);

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    int verified = verify(&p.piA, &p.piB2, &p.piC, p.uwProof, vk);

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
