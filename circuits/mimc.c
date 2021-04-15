#ifndef MIMC_C
#define MIMC_C

#define NROUNDS 91

void mimc7(element *h, element *x_in, element *k)
{
	char buff[2048];
	FILE *cnst;
	cnst = fopen("circuits/constants.txt", "r");

	element c[91];
	init_array(c, 91);

	for (int i = 0; i < 91; i++)
	{
		fgets(buff, sizeof buff, cnst);
		input(&c[i], buff);
	}

	fclose(cnst);

	element t[NROUNDS];
	element r[NROUNDS];
	element f[NROUNDS*3];

	init_array(t, NROUNDS);
	init_array(r, NROUNDS);
	init_array(f, NROUNDS*3);

	int it = 0;

	for (int i = 0; i < NROUNDS; i++)
	{
		if (i == 0) addmul(&t[i], k, x_in, &one);
		else add3mul(&t[i], k, &r[i-1], &c[i], &one);

		mul(&f[it], &t[i], &t[i]);
		mul(&f[it+1], &f[it], &f[it]);
		mul(&f[it+2], &f[it+1], &f[it]);
		mul(&r[i], &f[it+2], &t[i]);

		it = it + 3;
	}

	addmul(h, k, &r[NROUNDS-1], &one);
}

void multi_hash(element h, element *x_in, int arraySize)
{
	element k[arraySize];
	init_array(k, arraySize);
	input(&k[0], "0");

	for (int i = 0; i < arraySize-1; i++)
	{
		element hf;
		init(&hf);
		mimc7(&hf, &x_in[i], &k[i]);
		add3mul(&k[i+1], &x_in[i], &k[i], &hf, &one);
	}

	element hf;
	init(&hf);
	mimc7(&hf, &x_in[arraySize-1], &k[arraySize-1]);
	add3mul(&h, &x_in[arraySize-1], &k[arraySize-1], &hf, &one);
}

#endif