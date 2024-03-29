double elapsedSort;
double elapsedBosCoster;
char *transcript;

void generate_random_scalar(mclBnFr *value)
{
    if (test_no_rand) mclBnFr_setInt(value, 123456);
    else mclBnFr_setByCSPRNG(value);
}

void log_state(int type)
{
    if (logs)
    {
        if (type == 0) printf(" \033[1;31m[ERROR]\033[0m\n");
        if (type == 1) printf(" \033[1;32m[OK]\033[0m\n");
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
    if (type == 1 && bench) printf("\033[1;32m[SUCCESS] :\033[0m %s", msg);
    if (type == 0 && bench) printf("\033[1;31m[FAIL] :\033[0m %s", msg);
}

void init_setup(void *circuit)
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

    uw = (mpz_t*) malloc((M) * sizeof(mpz_t));
    LRO_constants = (mpz_t*) malloc((lro_const_total) * sizeof(mpz_t));

    for (int i = 0; i < M; i++)
    {
       mpz_init2(uw[i], BITS);
    }
}

void init_prover(void *circuit, proving_key pk)
{
    init_setup(circuit);

    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    int n = mpz_get_ui(pk.Ne);

    AsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    BsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));
    CsFr = (mclBnFr*) malloc((n) * sizeof(mclBnFr));

    mpz_init(pPrime);
    mpz_set_str(pPrime, PRIMESTR, 10);
    if (bench) printf("  |--- FFT domain size : %d\n", n);

    rsigma = (mpz_t*) malloc((n) * sizeof(mpz_t)); 
    rsigmaInv = (mpz_t*) malloc((n) * sizeof(mpz_t)); 

    mpz_t randNum;
    mpz_init(randNum);
    mpz_t factor;
    mpz_init_set_ui(factor, n);
    mpz_invert(factor, factor, pPrime);
    mpz_init(shift);

    mclBnFr rand;
    generate_random_scalar(&rand);
    fr_to_mpz(&randNum, &rand);
    mpz_set(shift, randNum);

    mpz_init2(rsigma[0], BITS);
    mpz_init2(rsigmaInv[0], BITS);
    mpz_set_ui(rsigma[0], 1);
    mpz_invert(rsigmaInv[0], rsigma[0], pPrime);

    mclBnFr frFactor;
    mpz_to_fr(&frFactor, &rsigmaInv[0]);
    mclBnG1_mul(&pk.xt1_rand[0], &pk.xt1[0], &frFactor);
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
        mclBnG1_mul(&pk.xt1_rand[i], &pk.xt1[i], &frFactorMulti);

        mpz_mul(rsigma[i], rsigma[i], factor);
        mpz_mod(rsigma[i], rsigma[i], pPrime);
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    log_success("ZPiE started successfully in ", 1);
    if (bench) printf("%fs\n\n", elapsed);
}

void bos_coster_bp(mclBnG1 *chunk, mclBnG1 *points, mclBnFr *scalars, int heapsize)
{
    mpz_t *exp[heapsize];
    mpz_t scalars_p[heapsize];
    mclBnG1 *points_p;
    points_p = (mclBnG1*) malloc((heapsize) * sizeof(mclBnG1));

    for (int i = 0; i < heapsize; i++)
    {
        fr_to_mpz(&scalars_p[i], &scalars[i]);
        points_p[i] = points[i];
        exp[i] = &scalars_p[i];
    }

    sort_list(exp, heapsize);
    while (mpz_cmp_ui(*exp[2], 0) != 0)
    {
        mpz_sub(*exp[0], *exp[0], *exp[2]); 
        mclBnG1_add(&points_p[exp[2]-scalars_p], &points_p[exp[0]-scalars_p], &points_p[exp[2]-scalars_p]);
        binarymaxheap(exp, 0, heapsize);
    }

    mclBnFr frFactor;
    mpz_to_fr(&frFactor, exp[0]);
    mclBnG1_mul(chunk, &points_p[exp[0]-scalars_p], &frFactor); 
}

static inline void mult_exp(mclBnG1 *chunk, mclBnG1 *points, mclBnFr *scalars, int heapsize)
{
    mclBnG1_mulVec(chunk, points, scalars, heapsize);
}

char *to_hex(const unsigned char *array, size_t length)
{
    char *outstr = malloc(2 * length + 1);
    if (!outstr) return outstr;

    char *p = outstr;
    for (size_t i = 0; i < length; i++) 
    {
        p += sprintf(p, "%02hhx", array[i]);
    }

    return outstr;
}

void transcript_hash(mclBnFr *hash)
{
    BYTE buff_hash_bytes[SHA256_BLOCK_SIZE];

    SHA256_CTX ctx;

    sha256_init(&ctx);
    //FIXME
    //sha256_update(&ctx, transcript, strlen(transcript));
    sha256_final(&ctx, buff_hash_bytes);

    char *buff_hash = to_hex(buff_hash_bytes, sizeof buff_hash_bytes);
    mclBnFr_setStr(hash, buff_hash, strlen(buff_hash)-1, 16);
}

static inline void transcript_add_Fr(mclBnFr *val)
{
    char buff[2048];
    mclBnFr_getStr(buff, sizeof(buff), val, 10);

    strcat(transcript, buff);
    strcat(transcript, "\n");
}

static inline void transcript_add_G1(mclBnG1 *val)
{
    char buff[2048];
    mclBnG1_getStr(buff, sizeof(buff), val, 10);

    strcat(transcript, buff);
    strcat(transcript, "\n");
}

void binarymaxheap(mpz_t *exp[], int i, int heapsize)
{
    int largest, left, right;
    mpz_t *temp;

    left = (2*i+1);
    right = ((2*i)+2);

    if (left >= heapsize) return;
    else
    {
        if (left < (heapsize) && (mpz_cmp(*exp[left], *exp[i]) > 0)) largest = left;
        else largest = i;
        if (right < (heapsize) && (mpz_cmp(*exp[right], *exp[largest]) > 0)) largest = right;
        if (largest != i)
        {
            temp = exp[i];
            exp[i] = exp[largest];
            exp[largest] = temp; 

            binarymaxheap(exp, largest, heapsize);
        }
    }
}

void sort_list(mpz_t *exp[], int heapsize)
{
    struct timespec begin, end;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    for (int j = heapsize/2; j >= 0; j--)
    {
        binarymaxheap(exp, j, heapsize);
    }   

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsedSort += (end.tv_sec - begin.tv_sec);
    elapsedSort += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
}

void bos_coster(mpz_t *exp[], int heapsize, int baseNum, proving_key *pk)
{
    sort_list(exp, heapsize);
    while (mpz_cmp_ui(*exp[2], 0) != 0)
    {
        struct timespec begin, end;
        clock_gettime(CLOCK_MONOTONIC, &begin);
        mpz_sub(*exp[0], *exp[0], *exp[2]); 

        if (baseNum) mclBnG1_add(&pk->xt1[exp[2]-wM], &pk->xt1[exp[0]-wM], &pk->xt1[exp[2]-wM]);
        else
        {
            mclBnG1_add(&pk->A1[exp[2]-uw], &pk->A1[exp[0]-uw], &pk->A1[exp[2]-uw]);
            mclBnG1_add(&pk->B1[exp[2]-uw], &pk->B1[exp[0]-uw], &pk->B1[exp[2]-uw]);
            mclBnG2_add(&pk->B2[exp[2]-uw], &pk->B2[exp[0]-uw], &pk->B2[exp[2]-uw]);
            mclBnG1_add(&pk->pk1[exp[2]-uw], &pk->pk1[exp[0]-uw], &pk->pk1[exp[2]-uw]);  
        }

        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsedBosCoster += (end.tv_sec - begin.tv_sec);
        elapsedBosCoster += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
        clock_gettime(CLOCK_MONOTONIC, &begin);
        binarymaxheap(exp, 0, heapsize);
        clock_gettime(CLOCK_MONOTONIC, &end);
        elapsedSort += (end.tv_sec - begin.tv_sec);
        elapsedSort += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
    }
}

int fr_cmp(mclBnFr *frFactor1, mclBnFr *frFactor2)
{
    mpz_t f1, f2;
    
    char buff[2048];
    mclBnFr_getStr(buff, sizeof(buff), frFactor1, 10);
    mpz_init_set_str(f1, buff, 10);
    mclBnFr_getStr(buff, sizeof(buff), frFactor2, 10);
    mpz_init_set_str(f2, buff, 10);

    return mpz_cmp(f1, f2);
}

void log_polynomial(mpz_t P[], int size, char letter[], int idx)
{
	if (logs)
    {
    	if (idx == -1) printf("%s(x) = ", letter);
    	else printf("%s%d(x) = ", letter, idx);

        for (int j = size-1; j >= 0; j--)
        {
            gmp_printf("%Zd", P[j]);
            if (j == 1)
            {
                printf("x");
            }
            else if (j>1)
            {
                printf("x^%d", j);
            }
            if (j > 0) printf(" + ");
        }
        printf("\n");
    }
}

int get_thread()
{
    #ifdef MULTI_SET
        return omp_get_thread_num();
    #else
        return 99;
    #endif
}

void mpz_to_fr(mclBnFr *frFactor, mpz_t *convert)
{
    char buff[2048];
    mpz_get_str(buff, 10, *convert);
    mclBnFr_setStr(frFactor, buff, strlen(buff), 10);
}

void fr_to_mpz(mpz_t *convert, mclBnFr *frFactor)
{
    char buff[2048];
    mclBnFr_getStr(buff, sizeof(buff), frFactor, 10);
    mpz_init(*convert);
    mpz_set_str(*convert, buff, 10);
}
