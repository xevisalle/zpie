mpz_t *u;

void generateqap(void *circuit, mpz_t *A, mpz_t *B, mpz_t *C, struct Trapdoor t)
{
    #pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mpz_init(A[i]);
        mpz_init(B[i]);
        mpz_init(C[i]);
    }

    L = (char **)malloc(N * sizeof(char*));
    R = (char **)malloc(N * sizeof(char*));
    O = (char **)malloc(N * sizeof(char*));

    for (int i = 0; i < N; i++)
    {
        L[i] = (char*) malloc(M * sizeof(char));
        R[i] = (char*) malloc(M * sizeof(char));
        O[i] = (char*) malloc(M * sizeof(char));
    }

    log_message("Computing R1CS...");
    init_circuit(circuit); 
    log_state(1);

    if (logs)
    {
        printf("A = \n");
        for (int i = 0; i < N; i++)
        {
            printf("( ");
            for (int j = 0; j < M; j++)
            {
                printf("%d ", L[i][j]);
            }

            printf(")\n");
        }

        printf("\nB = \n");
        for (int i = 0; i < N; i++)
        {
            printf("( ");
            for (int j = 0; j < M; j++)
            {
                printf("%d ", R[i][j]);
            }

            printf(")\n");
        }

        printf("\nC = \n");
        for (int i = 0; i < N; i++)
        {
            printf("( ");
            for (int j = 0; j < M; j++)
            {
                printf("%d ", O[i][j]);
            }

            printf(")\n");
        }
    }

    mpz_t factor, tX, T, uL;
    mpz_inits(factor, T, uL, NULL);
    mpz_init_set(tX, t.x);

    mpz_powm(T, tX, Ne, pPrime);
    mpz_sub_ui(T, T, 1); // Z^d - 1 = Zs

    mpz_invert(uL, Ne, pPrime);
    mpz_mul(uL, T, uL); // L1 = Zs / d
    mpz_mod(uL, uL, pPrime);

    u = (mpz_t*) malloc((N) * sizeof(mpz_t));

    for (int i = 0; i < N; i++)
    {
        mpz_init(u[i]);

        mpz_sub(factor, tX, wM[i]);
        mpz_invert(factor, factor, pPrime);
        mpz_mul(u[i], uL, factor);
        mpz_mod(u[i], u[i], pPrime);

        mpz_mul(uL, uL, wM[1]);
        mpz_mod(uL, uL, pPrime);
    }

    for (int j = N; j--;)
    {
        #pragma omp parallel for
        for (int i = 0; i < M; i++)
        {
            mpz_addmul_ui(A[i], u[j], L[j][i]);
            mpz_mod(A[i], A[i], pPrime);
            mpz_addmul_ui(B[i], u[j], R[j][i]);
            mpz_mod(B[i], B[i], pPrime);
            mpz_addmul_ui(C[i], u[j], O[j][i]);
            mpz_mod(C[i], C[i], pPrime);
        }
    }
    
    qapSize = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if(L[j][i]) qapSize += 3;
            if(R[j][i]) qapSize += 3;
            if(O[j][i]) qapSize += 3;
        }
    }
}