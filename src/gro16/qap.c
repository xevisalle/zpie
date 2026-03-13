
void generateqap(void* circuit, mclBnFr* A, mclBnFr* B, mclBnFr* C, struct Trapdoor t,
                 int* qap_size, int Ne)
{
#pragma omp parallel for
    for (int i = 0; i < M; i++)
    {
        mclBnFr_clear(&A[i]);
        mclBnFr_clear(&B[i]);
        mclBnFr_clear(&C[i]);
    }

    L = (int**) malloc(N * sizeof(int*));
    R = (int**) malloc(N * sizeof(int*));
    O = (int**) malloc(N * sizeof(int*));

    for (int i = 0; i < N; i++)
    {
        L[i] = (int*) malloc(M * sizeof(int));
        R[i] = (int*) malloc(M * sizeof(int));
        O[i] = (int*) malloc(M * sizeof(int));

        for (int j = 0; j < M; j++)
        {
            L[i][j] = 0;
            R[i][j] = 0;
            O[i][j] = 0;
        }
    }

    log_message("Computing R1CS...");

    cn = 0;
    lro_constants_n = 0;
    wn = nPublic + nConst;
    un = nConst;
    constant_n = 0;
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

    mclBnFr tX, T, uL, frNe, frOne, factor;
    tX = t.x;
    mclBnFr_setInt(&frNe, Ne);
    mclBnFr_setInt(&frOne, 1);

    mclBnFr_pow(&T, &tX, &frNe);
    mclBnFr_sub(&T, &T, &frOne); // Z^d - 1 = Zs

    mclBnFr_inv(&uL, &frNe);
    mclBnFr_mul(&uL, &T, &uL); // L1 = Zs / d

    mclBnFr* u = (mclBnFr*) malloc((N) * sizeof(mclBnFr));

    for (int i = 0; i < N; i++)
    {
        mclBnFr_sub(&factor, &tX, &wM[i]);
        mclBnFr_inv(&factor, &factor);
        mclBnFr_mul(&u[i], &uL, &factor);

        mclBnFr_mul(&uL, &uL, &wM[1]);
    }

    int l_it = lro_const_total - 2;
    int r_it = lro_const_total - 1;

    for (int j = N; j--;)
    {
        for (int i = M; i--;)
        {
            mclBnFr f;
            if (L[j][i] != INT_MAX)
            {
                mclBnFr frLji;
                mclBnFr_setInt(&frLji, L[j][i]);
                mclBnFr_mul(&f, &u[j], &frLji);
            }
            else
            {
                mclBnFr_mul(&f, &u[j], &LRO_constants[l_it]);
                l_it -= 2;
            }
            mclBnFr_add(&A[i], &A[i], &f);
            if (R[j][i] != INT_MAX)
            {
                mclBnFr frRji;
                mclBnFr_setInt(&frRji, R[j][i]);
                mclBnFr_mul(&f, &u[j], &frRji);
            }
            else
            {
                mclBnFr_mul(&f, &u[j], &LRO_constants[r_it]);
                r_it -= 2;
            }
            mclBnFr_add(&B[i], &B[i], &f);
            {
                mclBnFr frOji;
                mclBnFr_setInt(&frOji, O[j][i]);
                mclBnFr_mul(&f, &u[j], &frOji);
            }
            mclBnFr_add(&C[i], &C[i], &f);
        }
    }

    *qap_size = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (L[j][i] == 1)
                *qap_size += 3;
            else if (L[j][i] != 0)
                *qap_size += 4;
            if (R[j][i] == 1)
                *qap_size += 3;
            else if (R[j][i] != 0)
                *qap_size += 4;

            if (O[j][i])
                *qap_size += 3;
        }
    }

    free(u);
}