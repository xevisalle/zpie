
void addmul(element* oo, element* lo1, element* lo2, element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr_add(&uw[oo->index], &uw[lo1->index], &uw[lo2->index]);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &uw[ro->index]);
    }
    else
    {
        L[cn][lo1->index] = 1;
        L[cn][lo2->index] = 1;
        R[cn][ro->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void add3mul(element* oo, element* lo1, element* lo2, element* lo3, element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr_add(&uw[oo->index], &uw[lo1->index], &uw[lo2->index]);
        mclBnFr_add(&uw[oo->index], &uw[oo->index], &uw[lo3->index]);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &uw[ro->index]);
    }
    else
    {
        L[cn][lo1->index] = 1;
        L[cn][lo2->index] = 1;
        L[cn][lo3->index] = 1;
        R[cn][ro->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void addsmul(element* oo, int* size, element* los, element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        for (int i = 0; i < *size; i++)
        {
            mclBnFr_add(&uw[oo->index], &uw[oo->index], &uw[los[i].index]);
        }

        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &uw[ro->index]);
    }
    else
    {
        for (int i = 0; i < *size; i++)
        {
            L[cn][los[i].index] = 1;
        }

        R[cn][ro->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void add3muladd3(element* oo, element* lo1, element* lo2, element* lo3, element* ro1, element* ro2,
                 element* ro3)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr factor;
        mclBnFr_add(&uw[oo->index], &uw[lo1->index], &uw[lo2->index]);
        mclBnFr_add(&uw[oo->index], &uw[oo->index], &uw[lo3->index]);
        mclBnFr_add(&factor, &uw[ro1->index], &uw[ro2->index]);
        mclBnFr_add(&factor, &factor, &uw[ro3->index]);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &factor);
    }
    else
    {
        L[cn][lo1->index] = 1;
        L[cn][lo2->index] = 1;
        L[cn][lo3->index] = 1;
        R[cn][ro1->index] = 1;
        R[cn][ro2->index] = 1;
        R[cn][ro3->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void addmuladd(element* oo, element* lo1, element* lo2, element* ro1, element* ro2)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr factor;
        mclBnFr_add(&uw[oo->index], &uw[lo1->index], &uw[lo2->index]);
        mclBnFr_add(&factor, &uw[ro1->index], &uw[ro2->index]);
        mclBnFr_mul(&uw[oo->index], &factor, &uw[oo->index]);
    }
    else
    {
        L[cn][lo1->index] = 1;
        L[cn][lo2->index] = 1;
        R[cn][ro1->index] = 1;
        R[cn][ro2->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void mul(element* oo, element* lo, element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr_mul(&uw[oo->index], &uw[lo->index], &uw[ro->index]);
    }
    else
    {
        L[cn][lo->index] = 1;
        R[cn][ro->index] = 1;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void addmul_constants(element* oo, int* lc1, element* lo1, int* lc2, element* lo2, int* rc,
                      element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr factor, factor2;
        mclBnFr_setInt(&factor, *lc1);
        mclBnFr_mul(&factor, &uw[lo1->index], &factor);
        mclBnFr_setInt(&factor2, *lc2);
        mclBnFr_mul(&uw[oo->index], &uw[lo2->index], &factor2);
        mclBnFr_add(&factor, &factor, &uw[oo->index]);
        mclBnFr_setInt(&factor2, *rc);
        mclBnFr_mul(&uw[oo->index], &uw[ro->index], &factor2);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &factor);
    }
    else
    {
        L[cn][lo1->index] = *lc1;
        L[cn][lo2->index] = *lc2;
        R[cn][ro->index] = *rc;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void mul_constants(element* oo, int* lc, element* lo, int* rc, element* ro)
{
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr factor, factor2;
        mclBnFr_setInt(&factor, *lc);
        mclBnFr_mul(&factor, &uw[lo->index], &factor);
        mclBnFr_setInt(&factor2, *rc);
        mclBnFr_mul(&uw[oo->index], &uw[ro->index], &factor2);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &factor);
    }
    else
    {
        L[cn][lo->index] = *lc;
        R[cn][ro->index] = *rc;
        O[cn][oo->index] = 1;

        cn++;
    }
}

void mul_big_constants(element* oo, mclBnFr* lc, element* lo, mclBnFr* rc, element* ro)
{
    if (setParams)
    {
        lro_const_total += 2;
        N++;
    }
    else if (prover)
    {
        mclBnFr factor;
        mclBnFr_mul(&factor, &uw[lo->index], lc);
        mclBnFr_mul(&uw[oo->index], &uw[ro->index], rc);
        mclBnFr_mul(&uw[oo->index], &uw[oo->index], &factor);
    }
    else
    {
        L[cn][lo->index] = INT_MAX;
        R[cn][ro->index] = INT_MAX;
        O[cn][oo->index] = 1;

        cn++;
        mpz_t factor;
        fr_to_mpz(&factor, lc);
        mpz_set(LRO_constants[lro_constants_n], factor);
        fr_to_mpz(&factor, rc);
        mpz_set(LRO_constants[lro_constants_n + 1], factor);
        lro_constants_n += 2;
    }
}

void assert_equal(element* lo, element* ro)
{
    element factor1, factor2;
    init(&factor1);
    init(&factor2);

    mul(&factor1, ro, &oneNeg);
    
    if (setParams)
        N++;
    else if (prover)
    {
        mclBnFr_add(&uw[factor2.index], &uw[lo->index], &uw[factor1.index]);
    }
    else
    {
        L[cn][lo->index] = 1;
        L[cn][factor1.index] = 1;
        R[cn][one.index] = 1;

        cn++;
    }
}

void input(element* var, char* val)
{
    if (!setParams)
        mclBnFr_setStr(&uw[var->index], val, strlen(val), 10);
}

void init_constant(element* toAdd, char* val)
{
    if (setParams)
        M++;
    else
    {
        toAdd->index = constant_n;
        constant_n++;
        mclBnFr_setStr(&uw[toAdd->index], val, strlen(val), 10);
    }
    if (setParams)
        nConst++;
}

void init_public(element* toAdd)
{
    if (setParams)
        M++;
    else
    {
        toAdd->index = un;
        un++;
    }
    if (setParams)
        nPublic++;
}

void init_array(element* toAdd, int size)
{
    for (int i = 0; i < size; i++)
    {
        init(&toAdd[i]);
    }
}

void init(element* toAdd)
{
    if (setParams)
        M++;
    else
    {
        toAdd->index = wn;
        wn++;
    }
}

void init_circuit(void* circuit)
{
    init_constant(&one, "1");
    init_constant(&oneNeg, "-1");

    char buff[2048];
    FILE* cnst;
    cnst = fopen("gadgets/constants.txt", "r");

    for (int i = 0; i < 91; i++)
    {
        fgets(buff, sizeof buff, cnst);
        init_constant(&c_mimc[i], buff);
    }

    fclose(cnst);

    ((void (*)(void)) circuit)();
}
