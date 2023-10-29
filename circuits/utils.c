void add(element uOut, element vOut, element u1, element v1, element u2, element v2)
{
	element factor, factor1, factor2, factor3, factor4, factor5, factor6, factor7;
	init(&factor);
	init(&factor1);
	init(&factor2);
	init(&factor3);
	init(&factor4);
	init(&factor5);
	init(&factor6);
	init(&factor7);

	// uOut = (u1*v2 + v1*u2) / (1 + d*u1*u2*v1*v2)
	mul(&factor1, &u1, &v2);
	mul(&factor2, &v1, &u2);

	int d = 168696;
	int one_int = 1;

	mul_constants(&factor, &one_int, &factor1, &d, &factor2);

	mpz_t invFactor;
	mpz_init(invFactor);

	if(!setParams)
	{
		mpz_t f_check;
		mpz_init(f_check);
		mpz_add(f_check, uw[one.index], uw[factor.index]);
		mpz_invert(invFactor, f_check, pPrime);
	}

	char buff[2048];
  	mpz_get_str(buff, 10, invFactor);
  	input(&factor4, buff);

	addmul(&one, &factor, &one, &factor4); // verify x * 1/x = 1
  	addmul(&uOut, &factor1, &factor2, &factor4);

	// vOut = (v1*v2 - a*u1*u2) / (1 - d*u1*u2*v1*v2)
	mul(&factor5, &v1, &v2);
	
	int a = -168700;
	int one_neg = -1;

	mul_constants(&factor6, &a, &u1, &one_int, &u2);

	if(!setParams) 
	{
		mpz_t f_check;
		mpz_init(f_check);
		mpz_sub(f_check, uw[one.index], uw[factor.index]);
		mpz_invert(invFactor, f_check, pPrime);
	}

  	mpz_get_str(buff, 10, invFactor);
  	input(&factor7, buff);

	addmul_constants(&one, &one_int, &one, &one_neg, &factor, &one_int, &factor7); // verify x * 1/x = 1
  	addmul(&vOut, &factor5, &factor6, &factor7);
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

		mul_constants(&f4, &one_neg, &bits[i], &one_alone, &step1[i]);
		mul_constants(&f5, &one_neg, &bits[i], &one_alone, &step2[i]);

		if(i+1 != size)
		{
			add3mul(&step1[i+1], &f1, &f4, &step1[i], &one);
			add3mul(&step2[i+1], &f2, &f5, &step2[i], &one);
		}
		else
		{
			add3mul(&mulOut1, &f1, &f4, &step1[i], &one);
			add3mul(&mulOut2, &f2, &f5, &step2[i], &one);
		}
	}
}

void to_bits(element *bits, element val, int size)
{
	mpz_t t1, t2, t3, total;
	mpz_init(t1);
	mpz_init(t2);
	mpz_init(t3);
	mpz_init(total);

	mpz_set_str(t2, "1", 10);

	element b[size], fa;

	for (int i = 0; i < size; i++)
	{
		if(!setParams)
		{
			mpz_tdiv_q_2exp(t1, uw[val.index], i);
			mpz_and(t3, t1, t2);
		}

		char buff[2048];
		mpz_get_str(buff, 10, t3);
		input(&bits[i], buff);

		element power;
		mpz_ui_pow_ui(total, 2, i);
		mpz_get_str(buff, 10, total);
		init_constant(&power, buff);

		init(&fa);
		mul(&fa, &bits[i], &power);
		if (i == 0)
		{
			init(&b[i]);
			mul(&b[i], &fa, &one);
		}
		else
		{
			init(&b[i]);
			addmul(&b[i], &b[i-1], &fa, &one);
		}
	}

	assert_equal(&b[size-1], &val);
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