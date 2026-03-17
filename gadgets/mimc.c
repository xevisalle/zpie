#include <stdlib.h>

#ifndef MIMC_C
#define MIMC_C

#define NROUNDS 91

void mimc7(zpie_element* h, zpie_element* x_in, zpie_element* k)
{
    zpie_element r[NROUNDS];
    zpie_element f[NROUNDS * 3];

    zpie_init_array(r, NROUNDS);
    zpie_init_array(f, NROUNDS * 3);

    int it = 0;

    for (int i = 0; i < NROUNDS; i++)
    {
        if (i == 0)
            zpie_addmuladd(&f[it], k, x_in, k, x_in);
        else
            zpie_add3muladd3(&f[it], k, &r[i - 1], &c_mimc[i], k, &r[i - 1], &c_mimc[i]);

        zpie_mul(&f[it + 1], &f[it], &f[it]);
        zpie_mul(&f[it + 2], &f[it + 1], &f[it]);
        if (i == 0)
            zpie_addmul(&r[i], k, x_in, &f[it + 2]);
        else
            zpie_add3mul(&r[i], k, &r[i - 1], &c_mimc[i], &f[it + 2]);

        it = it + 3;
    }

    zpie_addmul(h, k, &r[NROUNDS - 1], &one);
}

void multi_hash(zpie_element h, zpie_element* x_in, int arraySize)
{
    zpie_element k[arraySize];
    zpie_init_array(k, arraySize);
    zpie_input(&k[0], "0");

    for (int i = 0; i < arraySize - 1; i++)
    {
        zpie_element hf;
        zpie_init(&hf);
        mimc7(&hf, &x_in[i], &k[i]);
        zpie_add3mul(&k[i + 1], &x_in[i], &k[i], &hf, &one);
    }

    zpie_element hf;
    zpie_init(&hf);
    mimc7(&hf, &x_in[arraySize - 1], &k[arraySize - 1]);
    zpie_add3mul(&h, &x_in[arraySize - 1], &k[arraySize - 1], &hf, &one);
}

#endif