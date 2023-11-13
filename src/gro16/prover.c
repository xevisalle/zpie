
void h_coefficients(proving_key *pk)
{
    int n = mpz_get_ui(pk->Ne);

    #pragma omp parallel for
    for (int j = 0; j < n; j++)
    {
        mclBnFr_clear(&AsFr[j]);
        mclBnFr_clear(&BsFr[j]);
        mclBnFr_clear(&CsFr[j]);
    }

    int l_it = 0;
    int r_it = 1;

    for (int j = 0; j < pk->qap_size; j+=3)
    {
        switch (pk->LRO[j])
        {
            case 1: mclBnFr_add(&AsFr[pk->LRO[j+1]], &AsFr[pk->LRO[j+1]], &uw[pk->LRO[j+2]]); break;
            case 2: mclBnFr_add(&BsFr[pk->LRO[j+1]], &BsFr[pk->LRO[j+1]], &uw[pk->LRO[j+2]]); break;
            case 3: mclBnFr_add(&CsFr[pk->LRO[j+1]], &CsFr[pk->LRO[j+1]], &uw[pk->LRO[j+2]]); break;
            case 10:
            {
                mclBnFr factorFr;
                if (pk->LRO[j+3] != INT_MAX)
                {
                    mclBnFr_setInt(&factorFr, pk->LRO[j+3]);
                    mclBnFr_mul(&factorFr, &uw[pk->LRO[j+2]], &factorFr);
                }
                else
                {
                    mclBnFr_mul(&factorFr, &uw[pk->LRO[j+2]], &pk->LRO_constants[l_it]);
                    l_it+=2;
                }
                mclBnFr_add(&AsFr[pk->LRO[j+1]], &AsFr[pk->LRO[j+1]], &factorFr); 
                j+=1;
                break;
            }
            case 20:
            {
                mclBnFr factorFr;
                if (pk->LRO[j+3] != INT_MAX)
                {
                    mclBnFr_setInt(&factorFr, pk->LRO[j+3]);
                    mclBnFr_mul(&factorFr, &uw[pk->LRO[j+2]], &factorFr);
                }
                else
                {
                    mclBnFr_mul(&factorFr, &uw[pk->LRO[j+2]], &pk->LRO_constants[r_it]);
                    r_it+=2;
                }
                mclBnFr_add(&BsFr[pk->LRO[j+1]], &BsFr[pk->LRO[j+1]], &factorFr);
                j+=1; 
                break;
            }
        }
    }

    #pragma omp parallel num_threads(3)
    {
        switch (get_thread())
        {
            case 0: ifft_t(n, pk->wM, AsFr); break;
            case 1: ifft_t(n, pk->wM, BsFr); break;
            case 2: ifft_t(n, pk->wM, CsFr); break;
            case 99:
                ifft_t(n, pk->wM, AsFr);
                ifft_t(n, pk->wM, BsFr);
                ifft_t(n, pk->wM, CsFr);
                break;
        }
    }

    #pragma omp parallel for
    for (int i = 0; i < n; i++)
    {       
        mclBnFr_mul(&AsFr[i], &AsFr[i], &BsFr[i]);
        mclBnFr_sub(&AsFr[i], &AsFr[i], &CsFr[i]);
    }

    ifft(n, pk->wM, AsFr);
}

void mul_exp(struct mulExpResult *result, mclBnFr *uwProof, proving_key *pk)
{
    int n = mpz_get_ui(pk->Ne);

    for (int i = nConst; i < (nPublic + nConst); i++)
    {
        uwProof[i-nConst] = uw[i];
    }

    #ifdef IS_MAC_OS
        int num_threads = 8;
    #else
        int num_threads = get_nprocs();
    #endif

    mclBnG1_mulVecMT(&result->uwA1, pk->A1, uw, M, num_threads);
    mclBnG1_mulVecMT(&result->uwB1, pk->B1, uw, M, num_threads);
    mclBnG2_mulVecMT(&result->uwB2, pk->B2, uw, M, num_threads);
    mclBnG1_mulVecMT(&result->uwC1, pk->pk1, uw + nPublic + nConst, M-(nPublic + nConst), num_threads);
    mclBnG1_mulVecMT(&result->htdelta, pk->xt1_rand, AsFr, n, num_threads);
}

void prove(int *circuit, mclBnG1 *piA, mclBnG2 *piB2, mclBnG1 *piC, mclBnFr *uwProof, proving_key *pk)
{
    prover = 1;

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    init_circuit(circuit);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    if (bench) printf("  |--- Circuit evaluation:  [%fs]\n", elapsed);

    prover = 0;

    clock_gettime(CLOCK_MONOTONIC, &begin);

    h_coefficients(pk);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    if (bench) printf("  |--- Compute h coefficients:  [%fs]\n", elapsed);

    clock_gettime(CLOCK_MONOTONIC, &begin);

    struct mulExpResult result;
    mul_exp(&result, uwProof, pk);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    if (bench) printf("  |--- G1, G2 multiexponentiations:  [%fs]\n", elapsed);
    log_message("Computing piA, piB1, piB2, piC, htdelta...");

    mclBnG1 piB1;
    mclBnFr r, s;

    generate_random_scalar(&r);
    generate_random_scalar(&s);

    // piA = s1.alpha + Auw + r * s1.delta;
    mclBnG1_mul(piA, &pk->delta1, &r);
    mclBnG1_add(piA, piA, &result.uwA1);
    mclBnG1_add(piA, piA, &pk->alpha1);
    // piB1 = s1.beta + B1uw + s * s1.delta;
    mclBnG1_mul(&piB1, &pk->delta1, &s);
    mclBnG1_add(&piB1, &piB1, &result.uwB1);
    mclBnG1_add(&piB1, &piB1, &pk->beta1);
    // piB2 = s2.beta + B2uw + s * s2.delta;
    mclBnG2_mul(piB2, &pk->delta2, &s);
    mclBnG2_add(piB2, piB2, &result.uwB2);
    mclBnG2_add(piB2, piB2, &pk->beta2);

    mclBnG1 factorG1;

    // piC = Cw + htdelta + piA*s + piB*r - r*s*s1.delta
    mclBnG1_mul(&factorG1, &pk->delta1, &r);
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