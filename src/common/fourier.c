void fft(size_t arr, mclBnFr domain[], mclBnFr *o)
{
    if (arr == 1) return;
    else
    {
        size_t arrNew = arr/2;

        mclBnFr oddDom[arrNew];
        mclBnFr oddVals[arrNew];
        mclBnFr evenVals[arrNew];

        for (int i = 0; i < arrNew; i++)
        {
            oddDom[i] = domain[2*i];
            oddVals[i] = o[2*i];
            evenVals[i] = o[(2*i)+1];
        }

        fft(arrNew, oddDom, oddVals);
        fft(arrNew, oddDom, evenVals);

        for (int i = 0; i < arrNew; i++)
        {
            mclBnFr_mul(&oddDom[i], &evenVals[i], &domain[i]);
            mclBnFr_add(&o[i], &oddVals[i], &oddDom[i]);
            mclBnFr_sub(&o[i+arrNew], &oddVals[i], &oddDom[i]);
        }
    }
}

void ifft(size_t arr, mclBnFr domain[], mclBnFr *o, mpz_t *Ne)
{
    fft(arr, domain, o);

    mpz_t factor, factor2;
    mclBnFr out[arr];
    mclBnFr frFactor;

    for (int i = 0; i < arr; i++)
    {
        out[i] = o[i];
    }

    mpz_init(factor2);
    mpz_powm(factor2, shift, *Ne, pPrime);
    mpz_sub_ui(factor2, factor2, 1);
    mpz_invert(factor2, factor2, pPrime);

    mpz_init_set_ui(factor, arr);
    mpz_invert(factor, factor, pPrime);
    mpz_mul(factor, factor, factor2);
    mpz_mod(factor, factor, pPrime);

    mpz_to_fr(&frFactor, &factor);
    mclBnFr_mul(&o[0], &out[0], &frFactor);

    for (int i = 1; i < arr; i++)
    {
        mclBnFr_mul(&o[i], &out[arr-i], &frFactor);
    }
}

void ifft_t(size_t arr, mclBnFr domain[], mclBnFr *o)
{
    fft(arr, domain, o);

    mclBnFr out[arr];

    for (int i = 0; i < arr; i++)
    {
        out[i] = o[i];
    }

    mclBnFr frFactor;
    mpz_to_fr(&frFactor, &rsigma[0]);
    mclBnFr_mul(&o[0], &out[0], &frFactor);

    for (int i = 1; i < arr; i++)
    {
        mpz_to_fr(&frFactor, &rsigma[i]);
        mclBnFr_mul(&o[i], &out[arr-i], &frFactor);
    }

    fft(arr, domain, o);
}