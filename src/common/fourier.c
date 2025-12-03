void fft(size_t arr, mclBnFr domain[], mclBnFr* o)
{
    if (arr == 1)
        return;
    else
    {
        size_t arrNew = arr / 2;

        mclBnFr oddDom[arrNew];
        mclBnFr oddVals[arrNew];
        mclBnFr evenVals[arrNew];

        for (int i = 0; i < arrNew; i++)
        {
            oddDom[i] = domain[2 * i];
            oddVals[i] = o[2 * i];
            evenVals[i] = o[(2 * i) + 1];
        }

        fft(arrNew, oddDom, oddVals);
        fft(arrNew, oddDom, evenVals);

        for (int i = 0; i < arrNew; i++)
        {
            mclBnFr_mul(&oddDom[i], &evenVals[i], &domain[i]);
            mclBnFr_add(&o[i], &oddVals[i], &oddDom[i]);
            mclBnFr_sub(&o[i + arrNew], &oddVals[i], &oddDom[i]);
        }
    }
}

void ifft(size_t arr, mclBnFr domain[], mclBnFr* o)
{
    fft(arr, domain, o);

    mclBnFr out[arr];
    mclBnFr frFactor;

    for (int i = 0; i < arr; i++)
    {
        out[i] = o[i];
    }

    mclBnFr_setInt(&frFactor, arr);
    mclBnFr_inv(&frFactor, &frFactor);
    mclBnFr_mul(&frFactor, &frFactor, &shift_fft);

    mclBnFr_mul(&o[0], &out[0], &frFactor);

    for (int i = 1; i < arr; i++)
    {
        mclBnFr_mul(&o[i], &out[arr - i], &frFactor);
    }
}

void ifft_t(size_t arr, mclBnFr domain[], mclBnFr* o)
{
    fft(arr, domain, o);

    mclBnFr out[arr];

    for (int i = 0; i < arr; i++)
    {
        out[i] = o[i];
    }

    mclBnFr_mul(&o[0], &out[0], &rsigma[0]);

    for (int i = 1; i < arr; i++)
    {
        mclBnFr_mul(&o[i], &out[arr - i], &rsigma[i]);
    }

    fft(arr, domain, o);
}