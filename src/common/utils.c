double elapsedSort;
double elapsedBosCoster;
char* transcript;

void generate_random_scalar(mclBnFr* value)
{
    if (test_no_rand)
        mclBnFr_setInt(value, 123456);
    else
        mclBnFr_setByCSPRNG(value);
}

void log_state(int type)
{
    if (logs)
    {
        if (type == 0)
            printf(" \033[1;31m[ERROR]\033[0m\n");
        if (type == 1)
            printf(" \033[1;32m[OK]\033[0m\n");
    }
}

void log_message(char msg[])
{
    if (logs)
    {
        printf("\033[1;34m[log] :\033[0m %s", msg);
    }
}

void log_success(char msg[], int type)
{
    if (type == 1 && bench)
        printf("\033[1;32m[SUCCESS] :\033[0m %s", msg);
    if (type == 0 && bench)
        printf("\033[1;31m[FAIL] :\033[0m %s", msg);
}

void init_setup(void* circuit)
{
    M = 0;
    N = 0;
    nPublic = 0;
    nConst = 0;
    lro_const_total = 0;

    mclBn_init(USEDCURVE, MCLBN_COMPILED_TIME_VAR);

    setParams = 1;
    init_circuit(circuit);
    setParams = 0;

    free(uw);
    free(LRO_constants);
    uw = (mclBnFr*) malloc((M) * sizeof(mclBnFr));
    LRO_constants = (mclBnFr*) malloc((lro_const_total) * sizeof(mclBnFr));

    for (int i = 0; i < M; i++)
    {
        mclBnFr_clear(&uw[i]);
    }
}

void init_prover(void* circuit, proving_key* pk)
{
    init_setup(circuit);

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    int n = pk->Ne;

    free(AsFr);
    free(BsFr);
    free(CsFr);
    AsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    BsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    CsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));

    if (bench)
        printf("  |--- FFT domain size : %d\n", n);

    free(rsigma);
    free(rsigmaInv);
    rsigma = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    rsigmaInv = (mclBnFr*) malloc((n) * sizeof(mclBnFr));

    mclBnFr rand, frNe, frOne;
    generate_random_scalar(&rand);
    mclBnFr_setInt(&frNe, n);
    mclBnFr_setInt(&frOne, 1);

    // shift = rand^Ne - 1, then invert
    mclBnFr shift;
    mclBnFr_pow(&shift, &rand, &frNe);
    mclBnFr_sub(&shift, &shift, &frOne);
    mclBnFr_inv(&shift, &shift);

    shift_fft = shift;

    mclBnFr_setInt(&rsigma[0], 1);
    mclBnFr_inv(&rsigmaInv[0], &rsigma[0]);

    mclBnG1_mul(&pk->xt1_rand[0], &pk->xt1[0], &rsigmaInv[0]);

    mclBnFr n_inverted;
    mclBnFr_setInt(&n_inverted, n);
    mclBnFr_inv(&n_inverted, &n_inverted);

    mclBnFr_mul(&rsigma[0], &rsigma[0], &n_inverted);

    mclBnFr one;
    mclBnFr_setInt(&one, 1);
    mclBnFr_mul(&rsigma[1], &rand, &one);

    for (int i = 1; i < n; i++)
    {
        if (i < n - 1)
            mclBnFr_mul(&rsigma[i + 1], &rsigma[i], &rand);

        mclBnFr_inv(&rsigmaInv[i], &rsigma[i]);
        mclBnG1_mul(&pk->xt1_rand[i], &pk->xt1[i], &rsigmaInv[i]);

        mclBnFr_mul(&rsigma[i], &rsigma[i], &n_inverted);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("ZPiE started successfully in ", 1);
    if (bench)
        printf("%fs\n\n", elapsed);
}

static inline void mult_exp(mclBnG1* chunk, mclBnG1* points, mclBnFr* scalars, int heapsize)
{
    mclBnG1_mulVec(chunk, points, scalars, heapsize);
}

char* to_hex(const unsigned char* array, size_t length)
{
    char* outstr = malloc(2 * length + 1);
    if (!outstr)
        return outstr;

    char* p = outstr;
    for (size_t i = 0; i < length; i++)
    {
        p += sprintf(p, "%02hhx", array[i]);
    }

    return outstr;
}

int get_thread()
{
#ifdef MULTI_SET
    return omp_get_thread_num();
#else
    return 99;
#endif
}
