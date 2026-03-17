
void zpie_addmul(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* ro)
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

void zpie_add3mul(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* lo3, zpie_element* ro)
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

void zpie_addsmul(zpie_element* oo, int* size, zpie_element* los, zpie_element* ro)
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

void zpie_add3muladd3(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* lo3, zpie_element* ro1, zpie_element* ro2,
                 zpie_element* ro3)
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

void zpie_addmuladd(zpie_element* oo, zpie_element* lo1, zpie_element* lo2, zpie_element* ro1, zpie_element* ro2)
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

void zpie_mul(zpie_element* oo, zpie_element* lo, zpie_element* ro)
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

void zpie_addmul_constants(zpie_element* oo, int* lc1, zpie_element* lo1, int* lc2, zpie_element* lo2, int* rc,
                      zpie_element* ro)
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

void zpie_mul_constants(zpie_element* oo, int* lc, zpie_element* lo, int* rc, zpie_element* ro)
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

void zpie_mul_big_constants(zpie_element* oo, mclBnFr* lc, zpie_element* lo, mclBnFr* rc, zpie_element* ro)
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
        LRO_constants[lro_constants_n] = *lc;
        LRO_constants[lro_constants_n + 1] = *rc;
        lro_constants_n += 2;
    }
}

void zpie_assert_equal(zpie_element* lo, zpie_element* ro)
{
    zpie_element factor1, factor2;
    zpie_init(&factor1);
    zpie_init(&factor2);

    zpie_mul(&factor1, ro, &oneNeg);

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

void zpie_input(zpie_element* var, char* val)
{
    if (!setParams)
        mclBnFr_setStr(&uw[var->index], val, strlen(val), 10);
}

void init_constant(zpie_element* toAdd, char* val)
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

void zpie_init_public(zpie_element* toAdd)
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

void zpie_init_array(zpie_element* toAdd, int size)
{
    for (int i = 0; i < size; i++)
    {
        zpie_init(&toAdd[i]);
    }
}

void zpie_init(zpie_element* toAdd)
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
        if (!fgets(buff, sizeof buff, cnst)) break;
        init_constant(&c_mimc[i], buff);
    }

    fclose(cnst);

    ((void (*)(void)) circuit)();
}
