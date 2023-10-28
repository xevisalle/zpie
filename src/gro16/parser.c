
void element_log(char *text, element *oo)
{
	if(!setParams) gmp_printf("%s%Zd\n", text, uw[oo->index]);
}

void addmul(element *oo, element *lo1, element *lo2, element *ro)
{
	if (setParams) N++;
	else if (prover)
	{
		mpz_add(uw[oo->index], uw[lo1->index], uw[lo2->index]);
		mpz_mul(uw[oo->index], uw[oo->index], uw[ro->index]);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
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

void addmul0(element *oo, element *lo1, element *lo2, element *ro)
{
	if (setParams) N++;
	else if (prover)
	{
		mpz_add(uw[oo->index], uw[lo1->index], uw[lo2->index]);
		mpz_mul(uw[oo->index], uw[oo->index], uw[ro->index]);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
	}
	else
	{
		L[cn][lo1->index] = 1;
		L[cn][lo2->index] = 1;
		R[cn][ro->index] = 1;

		cn++;
	}	
}

void add3mul(element *oo, element *lo1, element *lo2, element *lo3, element *ro)
{
	if (setParams) N++;
	else if (prover)
	{
		mpz_add(uw[oo->index], uw[lo1->index], uw[lo2->index]);
		mpz_add(uw[oo->index], uw[oo->index], uw[lo3->index]);
		mpz_mul(uw[oo->index], uw[oo->index], uw[ro->index]);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
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

void add3muladd3(element *oo, element *lo1, element *lo2, element *lo3, element *ro1, element *ro2, element *ro3)
{
	if (setParams) N++;
	else if (prover)
	{
		mpz_t factor;
		mpz_init(factor);
		mpz_add(uw[oo->index], uw[lo1->index], uw[lo2->index]);
		mpz_add(uw[oo->index], uw[oo->index], uw[lo3->index]);
		mpz_add(factor, uw[ro1->index], uw[ro2->index]);
		mpz_add(factor, factor, uw[ro3->index]);
		mpz_mul(uw[oo->index], uw[oo->index], factor);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
		mpz_clear(factor);
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

void addmuladd(element *oo, element *lo1, element *lo2, element *ro1, element *ro2)
{
	if (setParams) N++;
	else if (prover)
	{
		mpz_t factor;
		mpz_init(factor);
		mpz_add(uw[oo->index], uw[lo1->index], uw[lo2->index]);
		mpz_add(factor, uw[ro1->index], uw[ro2->index]);
		mpz_mul(uw[oo->index], factor, uw[oo->index]);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
		mpz_clear(factor);
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

void mul(element *oo, element *lo, element *ro)
{	
	if (setParams) N++;
	else if (prover)
	{
		mpz_mul(uw[oo->index], uw[lo->index], uw[ro->index]);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
	}
	else
	{
		L[cn][lo->index] = 1;
		R[cn][ro->index] = 1;
		O[cn][oo->index] = 1;

		cn++;
	}
}

void addmul_constants(element *oo, int *lc1, element *lo1, int *lc2, element *lo2, int *rc, element *ro)
{	
	if (setParams) N++;
	else if (prover)
	{
		mpz_t factor;
		mpz_init(factor);
		mpz_mul_si(factor, uw[lo1->index], *lc1);
		mpz_mul_si(uw[oo->index], uw[lo2->index], *lc2);
		mpz_add(factor, factor, uw[oo->index]);
		mpz_mul_si(uw[oo->index], uw[ro->index], *rc);
		mpz_mul(uw[oo->index], uw[oo->index], factor);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
		mpz_clear(factor);
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

void mul_constants(element *oo, int *lc, element *lo, int *rc, element *ro)
{	
	if (setParams) N++;
	else if (prover)
	{
		mpz_t factor;
		mpz_init(factor);
		mpz_mul_si(factor, uw[lo->index], *lc);
		mpz_mul_si(uw[oo->index], uw[ro->index], *rc);
		mpz_mul(uw[oo->index], uw[oo->index], factor);
		mpz_mod(uw[oo->index], uw[oo->index], pPrime);
		mpz_clear(factor);
	}
	else
	{
		L[cn][lo->index] = *lc;
		R[cn][ro->index] = *rc;
		O[cn][oo->index] = 1;

		cn++;
	}
}

void assert_equal(element *lo, element *ro)
{
	element factor1, factor2;
	init(&factor1);
	init(&factor2);

	mul(&factor1, ro, &oneNeg);
	addmul0(&factor2, lo, &factor1, &one);
}

void input(element *var, char *val)
{
	if (!setParams) mpz_set_str(uw[var->index], val, 10);
}

void init_public(element *toAdd)
{
	if (setParams) M++;
	else
	{
		toAdd->index = un;
		un++;
	}
	if (setParams) nPublic++;
}

void init_array(element *toAdd, int size)
{
	for (int i = 0; i < size; i++)
	{
		init(&toAdd[i]);
	}
}

void init(element *toAdd)
{
	if (setParams) M++;
	else
	{
		toAdd->index = wn;
		wn++;
	}
}

void init_circuit(void *circuit)
{
	init_public(&one);
	init_public(&oneNeg);

	input(&one, "1");
	input(&oneNeg, "-1");

	((void(*)(void))circuit)();
}

void test_full_api()
{
	element e_mul, e_addmul, e_add3mul, e_addmuladd;
    init(&e_mul);
	init(&e_addmul);
	init(&e_add3mul);
	init(&e_addmuladd);

    element a, b;
    init(&a);
    init(&b);

    input(&a, "5");
    input(&b, "10");

    mul(&e_mul, &a, &b);
	addmul(&e_addmul, &a, &b, &b);
	add3mul(&e_add3mul, &a, &a, &a, &b);
	addmuladd(&e_addmuladd, &a, &a, &b, &b);
}

void test_constraint_system(void)
{
	uw = (mpz_t*) malloc((8) * sizeof(mpz_t));
	wn = nPublic;
	un = 0;

    for (int i = 0; i < 8; i++)
    {
        mpz_init2(uw[i], BITS);
    }

	prover = 1;
	init_circuit(&test_full_api);
	prover = 0;

	CU_ASSERT(mpz_cmp_ui(uw[2], 50) == 0);
	CU_ASSERT(mpz_cmp_ui(uw[3], 150) == 0);
	CU_ASSERT(mpz_cmp_ui(uw[4], 150) == 0);
	CU_ASSERT(mpz_cmp_ui(uw[5], 200) == 0);
}