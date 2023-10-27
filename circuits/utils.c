void add(element uOut, element vOut, element u1, element v1, element u2, element v2)
{
	element a, d, dNeg;
	init(&a);
	init(&d);
	init(&dNeg);
	input(&a, "-168700");
	input(&d, "168696");
	input(&dNeg, "-168696");

	element factor, factor1, factor2, factor3, factor4, factor5, factor6, factor7, factor8, factor9, factor10;
	init(&factor);
	init(&factor1);
	init(&factor2);
	init(&factor3);
	init(&factor4);
	init(&factor5);
	init(&factor6);
	init(&factor7);
	init(&factor8);
	init(&factor9);
	init(&factor10);

	// uOut = (u1*v2 + v1*u2) / (1 + d*u1*u2*v1*v2)
	mul(&factor1, &u1, &v2);

	mul(&factor2, &v1, &u2);

	mul(&factor3, &factor1, &factor2);
	mul(&factor, &factor3, &d);

	addmul(&factor4, &factor, &one, &one);

	mpz_t invFactor;
	mpz_init(invFactor);
	if(!setParams) mpz_invert(invFactor, uw[factor4.index], pPrime);

	char buff[2048];
  	mpz_get_str(buff, 10, invFactor);
  	input(&factor5, buff);
  	addmul(&uOut, &factor1, &factor2, &factor5);

	// vOut = (v1*v2 - a*u1*u2) / (1 - d*u1*u2*v1*v2)
	mul(&factor6, &v1, &v2);
	mul(&factor7, &u1, &u2);
	mul(&factor8, &factor7, &a);

	element factorNeg;
	init(&factorNeg);

	mul(&factorNeg, &factor3, &dNeg);
	addmul(&factor9, &one, &factorNeg, &one);

	if(!setParams) mpz_invert(invFactor, uw[factor9.index], pPrime);
  	mpz_get_str(buff, 10, invFactor);
  	input(&factor10, buff);
  	addmul(&vOut, &factor6, &factor8, &factor10);
}

void mul_scalar(element mulOut1, element mulOut2, element A1, element A2, element *bits, int size)
{
	element accumulatedP1[size+1];
	element accumulatedP2[size+1];

	element step1[size+1];
	element step2[size+1];

	element doubledP1[size];
	element doubledP2[size];

	element oneNeg;
	init(&oneNeg);
	input(&oneNeg, "-1");

	init_array(doubledP1, size);
	init_array(doubledP2, size);

	init_array(accumulatedP1, size+1);
	init_array(accumulatedP2, size+1);
	init_array(step1, size+1);
	init_array(step2, size+1);

	input(&step1[0], "0");
	input(&step2[0], "1");

	int j;
	for (int i = 0; i < size; i++)
	{
		j = size - 1 - i;

		if (i == 0)
		{
			add(accumulatedP1[i+1], accumulatedP2[i+1], step1[i], step2[i], A1, A2);
			add(doubledP1[i], doubledP2[i], A1, A2, A1, A2);
		}
		else
		{
			add(accumulatedP1[i+1], accumulatedP2[i+1], step1[i], step2[i], doubledP1[i-1], doubledP2[i-1]);
			add(doubledP1[i], doubledP2[i], doubledP1[i-1], doubledP2[i-1], doubledP1[i-1], doubledP2[i-1]);
		}

		element f1, f2, f4, f5;
		init(&f1);
		init(&f2);
		init(&f4);
		init(&f5);

		mul(&f1, &accumulatedP1[i+1], &bits[i]);
		mul(&f2, &accumulatedP2[i+1], &bits[i]);

		int one_alone = 1;
		int one_neg = -1;

		addmul_constants(&f4, &one_neg, &oneNeg, &one_neg, &bits[i], &one_alone, &step1[i]);
		addmul_constants(&f5, &one_neg, &oneNeg, &one_neg, &bits[i], &one_alone, &step2[i]);

		addmul(&step1[i+1], &f1, &f4, &one);
		addmul(&step2[i+1], &f2, &f5, &one);
	}

	mul(&mulOut1, &step1[size], &one);
	mul(&mulOut2, &step2[size], &one);
}

void to_bits(element *bits, element val, int size)
{
	mpz_t t1, t2, t3, total;
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(total);

	element oneNeg;
	init(&oneNeg);
	input(&oneNeg, "-1");
	mpz_set_str(t2, "1", 10);

	for (int i = 0; i < size; i++)
	{
		if(!setParams)
		{
			mpz_tdiv_q_2exp(t1, uw[val.index], i);
			mpz_and(t3, t1, t2);
		}

		element b;
		init(&b);
		char buff[2048];
		mpz_get_str(buff, 10, t3);
		input(&bits[i], buff);

		addmul(&b, &bits[i], &oneNeg, &bits[i]);
		mpz_t pow;
		mpz_init(pow);
		mpz_ui_pow_ui(pow, 2, i);
		mpz_mul(t3, t3, pow);
		mpz_add(total, total, t3);
	}

	element check, checkCnst;
	init(&check);
	init(&checkCnst);
	char buff[2048];
	mpz_get_str(buff, 10, total);
	input(&check, buff);

	mul(&checkCnst, &check, &one);
}

typedef struct
{
    char *x; 
    char *y;
} point;

typedef struct
{
    point R; 
    char *S;
} eddsa_signature;