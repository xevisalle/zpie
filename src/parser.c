int prover;
int cn = 0;
int uwn = 0;

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

void input(element *var, char *val)
{
	if (!setParams) mpz_set_str(uw[var->index], val, 10);
}

void init_public(element *toAdd)
{
	init(toAdd);
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
		toAdd->index = uwn;
		uwn++;
	}
}

void init_circuit()
{
	init_public(&one);
	input(&one, "1");
	circuit();
}