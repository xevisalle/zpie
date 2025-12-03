#ifndef MIMC_C
#define MIMC_C

#define NROUNDS 91

void mimc7(element* h, element* x_in, element* k)
{
    element r[NROUNDS];
    element f[NROUNDS * 3];

    init_array(r, NROUNDS);
    init_array(f, NROUNDS * 3);

    int it = 0;

    for (int i = 0; i < NROUNDS; i++)
    {
        if (i == 0)
            addmuladd(&f[it], k, x_in, k, x_in);
        else
            add3muladd3(&f[it], k, &r[i - 1], &c_mimc[i], k, &r[i - 1], &c_mimc[i]);

        mul(&f[it + 1], &f[it], &f[it]);
        mul(&f[it + 2], &f[it + 1], &f[it]);
        if (i == 0)
            addmul(&r[i], k, x_in, &f[it + 2]);
        else
            add3mul(&r[i], k, &r[i - 1], &c_mimc[i], &f[it + 2]);

        it = it + 3;
    }

    addmul(h, k, &r[NROUNDS - 1], &one);
}

void multi_hash(element h, element* x_in, int arraySize)
{
    element k[arraySize];
    init_array(k, arraySize);
    input(&k[0], "0");

    for (int i = 0; i < arraySize - 1; i++)
    {
        element hf;
        init(&hf);
        mimc7(&hf, &x_in[i], &k[i]);
        add3mul(&k[i + 1], &x_in[i], &k[i], &hf, &one);
    }

    element hf;
    init(&hf);
    mimc7(&hf, &x_in[arraySize - 1], &k[arraySize - 1]);
    add3mul(&h, &x_in[arraySize - 1], &k[arraySize - 1], &hf, &one);
}

#endif