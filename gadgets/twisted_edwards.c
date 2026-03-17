void add_points(zpie_element uOut, zpie_element vOut, zpie_element u1, zpie_element v1, zpie_element u2, zpie_element v2)
{
    zpie_element factor, factor1, factor2, factor3, factor4, factor5, factor6, factor7;
    zpie_init(&factor);
    zpie_init(&factor1);
    zpie_init(&factor2);
    zpie_init(&factor3);
    zpie_init(&factor4);
    zpie_init(&factor5);
    zpie_init(&factor6);
    zpie_init(&factor7);

    // uOut = (u1*v2 + v1*u2) / (1 + d*u1*u2*v1*v2)
    zpie_mul(&factor1, &u1, &v2);
    zpie_mul(&factor2, &v1, &u2);

    int d = 168696;
    int one_int = 1;

    zpie_mul_constants(&factor, &one_int, &factor1, &d, &factor2);

    mclBnFr invFactor;

    if (!setParams)
    {
        mclBnFr f_check;
        mclBnFr_add(&f_check, &uw[one.index], &uw[factor.index]);
        mclBnFr_inv(&invFactor, &f_check);
    }

    char buff[2048];
    mclBnFr_getStr(buff, sizeof(buff), &invFactor, 10);
    zpie_input(&factor4, buff);

    zpie_addmul(&one, &factor, &one, &factor4); // verify x * 1/x = 1
    zpie_addmul(&uOut, &factor1, &factor2, &factor4);

    // vOut = (v1*v2 - a*u1*u2) / (1 - d*u1*u2*v1*v2)
    zpie_mul(&factor5, &v1, &v2);

    int a = -168700;
    int one_neg = -1;

    zpie_mul_constants(&factor6, &a, &u1, &one_int, &u2);

    if (!setParams)
    {
        mclBnFr f_check;
        mclBnFr_sub(&f_check, &uw[one.index], &uw[factor.index]);
        mclBnFr_inv(&invFactor, &f_check);
    }

    mclBnFr_getStr(buff, sizeof(buff), &invFactor, 10);
    zpie_input(&factor7, buff);

    zpie_addmul_constants(&one, &one_int, &one, &one_neg, &factor, &one_int,
                     &factor7); // verify x * 1/x = 1
    zpie_addmul(&vOut, &factor5, &factor6, &factor7);
}

void mul_scalar(zpie_element mulOut1, zpie_element mulOut2, zpie_element A1, zpie_element A2, zpie_element* bits, int size)
{
    zpie_element accumulatedP1[size + 1];
    zpie_element accumulatedP2[size + 1];

    zpie_element step1[size + 1];
    zpie_element step2[size + 1];

    zpie_element doubledP1[size];
    zpie_element doubledP2[size];

    zpie_init_array(doubledP1, size);
    zpie_init_array(doubledP2, size);

    zpie_init_array(accumulatedP1, size + 1);
    zpie_init_array(accumulatedP2, size + 1);
    zpie_init_array(step1, size + 1);
    zpie_init_array(step2, size + 1);

    zpie_input(&step1[0], "0");
    zpie_input(&step2[0], "1");

    int j;
    for (int i = 0; i < size; i++)
    {
        j = size - 1 - i;

        if (i == 0)
        {
            add_points(accumulatedP1[i + 1], accumulatedP2[i + 1], step1[i], step2[i], A1, A2);
            add_points(doubledP1[i], doubledP2[i], A1, A2, A1, A2);
        }
        else
        {
            add_points(accumulatedP1[i + 1], accumulatedP2[i + 1], step1[i], step2[i], doubledP1[i - 1],
                doubledP2[i - 1]);
            add_points(doubledP1[i], doubledP2[i], doubledP1[i - 1], doubledP2[i - 1], doubledP1[i - 1],
                doubledP2[i - 1]);
        }

        zpie_element f1, f2, f4, f5;
        zpie_init(&f1);
        zpie_init(&f2);
        zpie_init(&f4);
        zpie_init(&f5);

        zpie_mul(&f1, &accumulatedP1[i + 1], &bits[i]);
        zpie_mul(&f2, &accumulatedP2[i + 1], &bits[i]);

        int one_alone = 1;
        int one_neg = -1;

        zpie_mul_constants(&f4, &one_neg, &bits[i], &one_alone, &step1[i]);
        zpie_mul_constants(&f5, &one_neg, &bits[i], &one_alone, &step2[i]);

        if (i + 1 != size)
        {
            zpie_add3mul(&step1[i + 1], &f1, &f4, &step1[i], &one);
            zpie_add3mul(&step2[i + 1], &f2, &f5, &step2[i], &one);
        }
        else
        {
            zpie_add3mul(&mulOut1, &f1, &f4, &step1[i], &one);
            zpie_add3mul(&mulOut2, &f2, &f5, &step2[i], &one);
        }
    }
}

void to_bits(zpie_element* bits, zpie_element val, int size)
{
    unsigned char bytes[SIZE_FR];
    zpie_element b[size];

    for (int i = 0; i < size; i++)
    {
        int bit = 0;
        if (!setParams)
        {
            mclBnFr_getLittleEndian(bytes, SIZE_FR, &uw[val.index]);
            int byte_idx = i / 8;
            int bit_idx = i % 8;
            bit = (bytes[byte_idx] >> bit_idx) & 1;
        }

        char buff[2];
        buff[0] = '0' + bit;
        buff[1] = '\0';
        zpie_input(&bits[i], buff);

        mclBnFr one_mcl;
        mclBnFr_setInt(&one_mcl, 1);

        zpie_init(&b[i]);
        mclBnFr factor;
        mclBnFr_setInt(&factor, 1);
        // factor = 2^i
        for (int j = 0; j < i; j++)
        {
            mclBnFr two;
            mclBnFr_setInt(&two, 2);
            mclBnFr_mul(&factor, &factor, &two);
        }
        zpie_mul_big_constants(&b[i], &factor, &bits[i], &one_mcl, &one);
    }

    zpie_element fa;
    zpie_init(&fa);

    zpie_addsmul(&fa, &size, b, &one);
    zpie_assert_equal(&fa, &val);
}

typedef struct
{
    char* x;
    char* y;
} point;
