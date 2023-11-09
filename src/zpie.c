
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

    setup_keys keys;
    mpz_init_set(keys.pk.Ne, Ne);

    keys.pk.LRO_constants = (mpz_t*) malloc((lro_const_total) * sizeof(mpz_t));
    keys.pk.wMFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));

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
        mpz_to_fr(&keys.pk.wMFr[i], &wM[i]);
    }

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    setup(circuit, &t, &s1, &s2, &alphabetaT, &keys.pk.qap_size, &keys.pk.Ne);

    keys.pk.LRO = (int*) malloc((keys.pk.qap_size) * sizeof(int));

    for (int i = 0; i < lro_const_total; i++)
    {
        mpz_init_set(keys.pk.LRO_constants[i], LRO_constants[i]);
    }

    int it = 0;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(L[j][i] != 0) 
            {
                keys.pk.LRO[it+1] = j;
                keys.pk.LRO[it+2] = i;

                if (L[j][i] != 1)
                {
                    keys.pk.LRO[it] = 10;
                    keys.pk.LRO[it+3] = L[j][i];
                    it+=4;
                }
                else
                {
                    keys.pk.LRO[it] = 1;
                    it+=3;
                }
            }
            if(R[j][i] != 0) 
            {
                keys.pk.LRO[it+1] = j;
                keys.pk.LRO[it+2] = i;

                if (R[j][i] != 1)
                {
                    keys.pk.LRO[it] = 20;
                    keys.pk.LRO[it+3] = R[j][i];
                    it+=4;
                }
                else 
                {
                    keys.pk.LRO[it] = 2;
                    it+=3;
                }
            }
            if(O[j][i]) 
            {
                keys.pk.LRO[it] = 3;
                keys.pk.LRO[it+1] = j;
                keys.pk.LRO[it+2] = i;
                it+=3;
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

    keys.vk.constants = (mpz_t*) malloc((nConst) * sizeof(mpz_t));

    for (int i = 0; i < (nConst); i++)
    {
        mpz_init(keys.vk.constants[i]);
        mpz_set(keys.vk.constants[i], uw[i]);
    }

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        keys.vk.vk1[i] = s1.vk[i];
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("Setup generated successfully in", 1);
    if (bench) printf(" %fs\n", elapsed);

    return keys;
}

void serialize_pk(proving_key *pk)
{
    FILE *fpk;
    fpk = fopen("data/provingkey.params", "w");

    int n = mpz_get_ui(pk->Ne);

    int buff_pk_size = SIZE_FR * n + SIZE_G2 * (2 + M) + SIZE_G1 * (M - (nPublic + nConst) + 3 + n + 2 * M);
    char buff_pk[buff_pk_size];

    mpz_out_raw(fpk, pk->Ne);

    mpz_t factor;
    mpz_init(factor);
    mpz_set_si(factor, pk->qap_size);
    mpz_out_raw(fpk, factor);

    for (int i = 0; i < pk->qap_size; i++)
    {
        mpz_set_si(factor, pk->LRO[i]);
        mpz_out_raw(fpk, factor);
    }

    for (int i = 0; i < lro_const_total; i++)
    {
        mpz_out_raw(fpk, pk->LRO_constants[i]);
    }

    int size = 0;

    for (int i = 0; i < n; i++)
    {
        size += mclBnFr_serialize(buff_pk + size, SIZE_FR, &pk->wMFr[i]);
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

    for (int i = 0; i < M-(nPublic + nConst); i++)
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

void serialize_vk(verifying_key *vk)
{
    FILE *fvk;
    fvk = fopen("data/verifyingkey.params", "w");

    int buff_vk_size = SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst);
    char buff_vk[buff_vk_size];

    for (int i = 0; i < nConst; i++)
    {
        mpz_out_raw(fvk, vk->constants[i]);
    }

    int size = 0;

    size += mclBnGT_serialize(buff_vk, SIZE_GT, &vk->alphabetaT);
    size += mclBnG2_serialize(buff_vk + size, SIZE_G2, &vk->gamma2);
    size += mclBnG2_serialize(buff_vk + size, SIZE_G2, &vk->delta2);

    for (int i = 0; i < (nPublic + nConst); i++)
    {
        size += mclBnG1_serialize(buff_vk + size, SIZE_G1, &vk->vk1[i]);
    }

    fwrite(buff_vk, 1, size, fvk);
    fclose(fvk);
}

void store_setup(setup_keys *keys)
{
    struct stat st = {0};
    if (stat("data", &st) == -1) mkdir("data", 0700);

    serialize_pk(&keys->pk);
    serialize_vk(&keys->vk);
}

setup_keys read_setup(void *circuit)
{
    init_setup(circuit);

    FILE *fpk, *fvk;
    
    fpk = fopen("data/provingkey.params", "r");
    fvk = fopen("data/verifyingkey.params", "r");

    setup_keys keys;

    mpz_init(keys.pk.Ne);
    mpz_inp_raw(keys.pk.Ne, fpk);

    int n = mpz_get_ui(keys.pk.Ne);

    int buff_pk_size = SIZE_FR * n + SIZE_G2 * (2 + M) + SIZE_G1 * (M - (nPublic + nConst) + 3 + n + 2 * M);
    char buff_pk[buff_pk_size];
    
    keys.pk.wMFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    keys.vk.vk1 = (mclBnG1*) malloc(((nPublic + nConst)) * sizeof(mclBnG1));
    keys.vk.constants = (mpz_t*) malloc(((nConst)) * sizeof(mpz_t));

    keys.pk.xt1 = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.xt1_rand = (mclBnG1*) malloc((n) * sizeof(mclBnG1));
    keys.pk.A1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.B1 = (mclBnG1*) malloc((M) * sizeof(mclBnG1));
    keys.pk.pk1 = (mclBnG1*) malloc((M-(nPublic + nConst)) * sizeof(mclBnG1));
    keys.pk.B2 = (mclBnG2*) malloc((M) * sizeof(mclBnG2));
    keys.pk.LRO_constants = (mpz_t*) malloc((lro_const_total) * sizeof(mpz_t));

    mpz_t factor;
    mpz_init(factor);
    mpz_inp_raw(factor, fpk);

    keys.pk.qap_size = mpz_get_si(factor);
    keys.pk.LRO = (int*) malloc((keys.pk.qap_size) * sizeof(int)); 

    for (int i = 0; i < keys.pk.qap_size; i++)
    {
        mpz_inp_raw(factor, fpk);
        keys.pk.LRO[i] = mpz_get_si(factor);
    }

    for (int i = 0; i < lro_const_total; i++)
    {
        mpz_init(keys.pk.LRO_constants[i]);
        mpz_inp_raw(keys.pk.LRO_constants[i], fpk);
    }

    int size = 0;
    fread(buff_pk, 1, buff_pk_size, fpk);

    for (int i = 0; i < n; i++)
    {
        size += mclBnFr_deserialize(&keys.pk.wMFr[i], buff_pk + size, SIZE_FR);
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

    for (int i = 0; i < M-(nPublic + nConst); i++)
    {
        size += mclBnG1_deserialize(&keys.pk.pk1[i], buff_pk + size, SIZE_G1);
    }

    for (int i = 0; i < n; i++)
    {
        size += mclBnG1_deserialize(&keys.pk.xt1[i], buff_pk + size, SIZE_G1);
    }

    for (int i = 0; i < nConst; i++)
    {
        mpz_init(keys.vk.constants[i]);
        mpz_inp_raw(keys.vk.constants[i], fvk);
    }
    
    char buff_vk[SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst)];
    size = 0;

    fread(buff_vk, 1, SIZE_GT + SIZE_G2 * 2 + SIZE_G1 * (nPublic + nConst), fvk);
    size += mclBnGT_deserialize(&keys.vk.alphabetaT, buff_vk, SIZE_GT);
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

proof generate_proof(void *circuit, proving_key *pk)
{
    init_prover(circuit, pk);

    wn = nPublic + nConst;
    un = nConst;
    constant_n = 0;
    for (int i = 0; i < M; i++)
    {
        mpz_init(uw[i]);
    }

    int n = mpz_get_ui(pk->Ne);
    wM = (mpz_t*) malloc((n) * sizeof(mpz_t));
    
    proof p;

    p.uwProof = (mclBnFr*) malloc((nPublic) * sizeof(mclBnFr));

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

        mclBnFr_clear(&rsigma[i]);
        mclBnFr_clear(&rsigmaInv[i]);
    }

    for (int i = 0; i < M; i++)
    {
        mpz_clear(uw[i]);
    }

    return p;
}

void store_proof(proof *p)
{
    char buff[2048];
    FILE *fproof;
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
    FILE *fproof;
    fproof = fopen("data/proof.params", "r");

    p.uwProof = (mclBnFr*) malloc((nPublic) * sizeof(mclBnFr));

    int size = 0;

    fread(buff, 1, (SIZE_FR * nPublic) + SIZE_G1 + SIZE_G2 + SIZE_G1, fproof);

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

int verify_proof(void *circuit, proof *p, verifying_key *vk)
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
        if (bench) printf("%fs\n", elapsed);
    }
    else
    {
        log_success("Proof cannot be verified\n", 0);
    }

    return verified;
}
