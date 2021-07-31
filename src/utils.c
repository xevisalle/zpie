double elapsedSort;
double elapsedBosCoster;
char *transcript;

void bos_coster_bp(mclBnG1 *chunk, mclBnG1 *points, mclBnFr *scalars, int heapsize)
{
    mpz_t *exp[heapsize];
    mpz_t scalars_p[heapsize];
    mclBnG1 *points_p;
    points_p = (mclBnG1*) malloc((heapsize) * sizeof(mclBnG1));
    mclBnFr one;
    mclBnFr_setInt(&one, 1);

    for (int i = 0; i < heapsize; i++)
    {
        fr_to_mpz(&scalars_p[i], &scalars[i]);
        mclBnG1_mul(&points_p[i], &points[i], &one);
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
    #ifdef BOSCOSTER_MULEXP
        if ((heapsize >= 32) && (heapsize == Nb*Mc)) bos_coster_bp(chunk, points, scalars, heapsize);
        else mclBnG1_mulVec(chunk, points, scalars, heapsize);
    #elif MCL_MULEXP
        mclBnG1_mulVec(chunk, points, scalars, heapsize);
    #endif
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
    sha256_update(&ctx, transcript, strlen(transcript));
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

void bos_coster(mpz_t *exp[], int heapsize, int baseNum, struct ProvingKey *pk)
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
    #ifdef MULTI
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