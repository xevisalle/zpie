#define Nb 64
#define Mc 2

void bulletproof()
{
    // SETUP
    char buff1[2048], buff2[2048];
    mclBnG1 Gen;
    mclBnG1_setStr(&Gen, GGEN, strlen(GGEN), 10);

    mclBnFr rnd, *gamma;
    mclBnG1 Gb, Hb, Ub;
    mclBnG1 *G, *H;

    gamma = (mclBnFr*) malloc((Mc) * sizeof(mclBnFr));
    G = (mclBnG1*) malloc((Nb*Mc) * sizeof(mclBnG1));
    H = (mclBnG1*) malloc((Nb*Mc) * sizeof(mclBnG1));

    for (int i = 0; i < Mc; i++)
    {
        mclBnFr_setByCSPRNG(&gamma[i]);
    }
    
    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Gb, &Gen, &rnd);

    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Hb, &Gen, &rnd);

    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Ub, &Gen, &rnd);

    mclBnFr frFactor, frFactor2, frFactor3, frFactor4;

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_setInt(&frFactor, i);
        mclBnG1_mul(&G[i], &Gb, &frFactor);
        mclBnG1_mul(&H[i], &Hb, &frFactor);
    }

    //PROVER
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    mclBnFr *aL, *aR, *l, *r;
    mclBnFr *two_vec, v[Mc];

    aL = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    aR = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    l = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    r = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    two_vec = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));

    mclBnG1 g1Factor;
    mclBnG1 *V = (mclBnG1*) malloc((Mc) * sizeof(mclBnG1));

    mclBnFr_setInt(&two_vec[0], 1);
    mclBnFr_setInt(&two_vec[1], 2);

    for (int i = 2; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&two_vec[i], &two_vec[i-1], &two_vec[1]);
    }

    mclBnFr_setInt(&aL[3], 1);

    mclBnFr one;
    mclBnFr_setInt(&one, 1);
    mclBnFr alpha;
    mclBnFr_setByCSPRNG(&alpha);

    mclBnFr sL[Nb*Mc], sR[Nb*Mc], rho;
    mclBnFr_setByCSPRNG(&rho);

    mclBnG1 A;
    mclBnG1_mul(&A, &Hb, &alpha);
    mclBnG1 S;
    mclBnG1_mul(&S, &Hb, &rho);

    // compute inner products v[]
    for (int i = 0; i < Mc; i++)
    {
        mclBnFr_clear(&v[i]);
        for (int j = 0; j < Nb; j++)
        {
            mclBnFr_mul(&frFactor, &aL[j + (i * Nb)], &two_vec[j]);
            mclBnFr_add(&v[i], &v[i], &frFactor);
        }

        // compute commitment V
        mclBnG1_mul(&g1Factor, &Gb, &v[i]);
        mclBnG1_mul(&V[i], &Hb, &gamma[i]);
        mclBnG1_add(&V[i], &V[i], &g1Factor);
    }

    for (int i = 0; i < Nb*Mc; i++)
    {
        // compute aR
        mclBnFr_sub(&aR[i], &aL[i], &one);
        // compute commitment A
        mclBnG1_mul(&g1Factor, &G[i], &aL[i]);
        mclBnG1_add(&A, &A, &g1Factor);
        mclBnG1_mul(&g1Factor, &H[i], &aR[i]);
        mclBnG1_add(&A, &A, &g1Factor);
        // compute commitment S
        mclBnFr_setByCSPRNG(&sL[i]);
        mclBnFr_setByCSPRNG(&sR[i]);
        mclBnG1_mul(&g1Factor, &G[i], &sL[i]);
        mclBnG1_add(&S, &S, &g1Factor);
        mclBnG1_mul(&g1Factor, &H[i], &sR[i]);
        mclBnG1_add(&S, &S, &g1Factor);
    }

    mclBnFr y, z, z2;
    mclBnFr_setByCSPRNG(&y);
    mclBnFr_setByCSPRNG(&z);
    mclBnFr_mul(&z2, &z, &z);

    mclBnFr y_vec[Nb*Mc];
    mclBnFr_setInt(&y_vec[0], 1);
    mclBnFr_mul(&y_vec[1], &y_vec[0], &y);

    for (int i = 2; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&y_vec[i], &y_vec[i-1], &y);
    }

    mclBnFr t1, t2;
    mclBnFr_clear(&t1);
    mclBnFr_clear(&t2);
    mclBnFr_clear(&frFactor3);
    mclBnFr_mul(&frFactor3, &one, &z);
    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_add(&frFactor, &aR[i], &z);
        mclBnFr_mul(&frFactor, &y_vec[i], &frFactor);

        if (i % Nb == 0) mclBnFr_mul(&frFactor3, &frFactor3, &z);
        mclBnFr_mul(&frFactor2, &two_vec[i % Nb], &frFactor3);

        mclBnFr_add(&frFactor, &frFactor, &frFactor2);

        mclBnFr_mul(&frFactor, &sL[i], &frFactor);
        mclBnFr_add(&t1, &t1, &frFactor);

        mclBnFr_sub(&frFactor, &aL[i], &z);
        mclBnFr_mul(&frFactor2, &y_vec[i], &sR[i]);
        mclBnFr_mul(&frFactor, &frFactor, &frFactor2);
        mclBnFr_add(&t1, &t1, &frFactor);

        mclBnFr_mul(&frFactor, &y_vec[i], &sR[i]);
        mclBnFr_mul(&frFactor, &sL[i], &frFactor);
        mclBnFr_add(&t2, &t2, &frFactor);
    }

    mclBnFr tau1, tau2;
    mclBnG1 T1, T2;
    mclBnFr_setByCSPRNG(&tau1);
    mclBnFr_setByCSPRNG(&tau2);

    mclBnG1_mul(&g1Factor, &Gb, &t1);
    mclBnG1_mul(&T1, &Hb, &tau1);
    mclBnG1_add(&T1, &T1, &g1Factor);

    mclBnG1_mul(&g1Factor, &Gb, &t2);
    mclBnG1_mul(&T2, &Hb, &tau2);
    mclBnG1_add(&T2, &T2, &g1Factor);

    mclBnFr x, t_inner, tx, mu;
    mclBnFr_setByCSPRNG(&x);

    mclBnFr_clear(&frFactor2);
    mclBnFr_mul(&frFactor2, &one, &z);

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_sub(&l[i], &aL[i], &z);
        mclBnFr_mul(&frFactor, &sL[i], &x);
        mclBnFr_add(&l[i], &l[i], &frFactor);

        mclBnFr_mul(&frFactor, &sR[i], &x);
        mclBnFr_add(&frFactor, &frFactor, &z);
        mclBnFr_add(&frFactor, &frFactor, &aR[i]);
        mclBnFr_mul(&frFactor, &frFactor, &y_vec[i]);

        if (i % Nb == 0) mclBnFr_mul(&frFactor2, &frFactor2, &z);
        mclBnFr_mul(&r[i], &frFactor2, &two_vec[i % Nb]);
        mclBnFr_add(&r[i], &r[i], &frFactor);

        mclBnFr_mul(&frFactor, &l[i], &r[i]);
        mclBnFr_add(&t_inner, &t_inner, &frFactor);
    }

    mclBnFr_mul(&frFactor, &x, &x);
    mclBnFr_mul(&frFactor, &frFactor, &tau2);

    mclBnFr_mul(&tx, &tau1, &x);
    mclBnFr_add(&tx, &tx, &frFactor);

    mclBnFr_mul(&frFactor2, &one, &z);
    mclBnFr_clear(&frFactor3);
    for (int i = 1; i < Mc + 1; i++)
    {
        mclBnFr_mul(&frFactor2, &frFactor2, &z);
        mclBnFr_mul(&frFactor, &frFactor2, &gamma[i-1]);
        mclBnFr_add(&frFactor3, &frFactor3, &frFactor);
    }

    mclBnFr_add(&tx, &tx, &frFactor3);

    mclBnFr_mul(&mu, &rho, &x);
    mclBnFr_add(&mu, &mu, &alpha);

    mclBnG1 Hs[Nb*Mc], P;

    mclBnG1_mul(&P, &S, &x);
    mclBnG1_add(&P, &P, &A);

    mclBnFr_clear(&frFactor3);
    mclBnFr_mul(&frFactor3, &one, &z);

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_inv(&frFactor, &y_vec[i]);
        mclBnG1_mul(&Hs[i], &H[i], &frFactor);

        mclBnFr_mul(&frFactor, &z, &y_vec[i]);
        if (i % Nb == 0) mclBnFr_mul(&frFactor3, &frFactor3, &z);
        mclBnFr_mul(&frFactor2, &frFactor3, &two_vec[i % Nb]);
        mclBnFr_add(&frFactor, &frFactor, &frFactor2);

        mclBnG1_mul(&g1Factor, &Hs[i], &frFactor);
        mclBnG1_add(&P, &P, &g1Factor);

        mclBnFr_neg(&frFactor, &z);
        mclBnG1_mul(&g1Factor, &G[i], &frFactor);
        mclBnG1_add(&P, &P, &g1Factor);
    }

    mclBnG1 Pp;
    mclBnFr xp;
    mclBnFr_neg(&frFactor, &mu);
    mclBnG1_mul(&g1Factor, &Hb, &frFactor);
    mclBnG1_add(&Pp, &P, &g1Factor);

    // inner product
    mclBnFr_setByCSPRNG(&xp);
    mclBnFr_mul(&frFactor, &xp, &t_inner);
    mclBnG1_mul(&g1Factor, &Ub, &frFactor);
    mclBnG1_add(&Pp, &Pp, &g1Factor);

    mclBnG1_mul(&Ub, &Ub, &xp);

    int logN = log(Nb*Mc)/log(2);
    mclBnFr xp_vec[logN];
    mclBnG1 L_vec[logN], R_vec[logN];

    mclBnFr lp[Nb*Mc], rp[Nb*Mc];
    mclBnG1 Gp[Nb*Mc], Hp[Nb*Mc];

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&lp[i], &l[i], &one);
        mclBnFr_mul(&rp[i], &r[i], &one);
        mclBnG1_mul(&Gp[i], &G[i], &one);
        mclBnG1_mul(&Hp[i], &Hs[i], &one);
    }

    int np = Nb*Mc;
    mclBnFr cl, cr;
    for (int i = 0; i < logN; i++)
    {
        np = np / 2;

        mclBnFr_clear(&cl);
        mclBnFr_clear(&cr);

        for (int j = 0; j < np; j++)
        {
            mclBnFr_mul(&frFactor, &lp[j], &rp[j+np]);
            mclBnFr_add(&cl, &cl, &frFactor);
            mclBnFr_mul(&frFactor, &lp[j+np], &rp[j]);
            mclBnFr_add(&cr, &cr, &frFactor);
        }

        mclBnG1_mul(&L_vec[i], &Ub, &cl);
        mclBnG1_mul(&R_vec[i], &Ub, &cr);

        for (int j = 0; j < np; j++)
        {
            mclBnG1_mul(&g1Factor, &Gp[j+np], &lp[j]);
            mclBnG1_add(&L_vec[i], &L_vec[i], &g1Factor);
            mclBnG1_mul(&g1Factor, &Hp[j], &rp[j+np]);
            mclBnG1_add(&L_vec[i], &L_vec[i], &g1Factor);

            mclBnG1_mul(&g1Factor, &Gp[j], &lp[j+np]);
            mclBnG1_add(&R_vec[i], &R_vec[i], &g1Factor);
            mclBnG1_mul(&g1Factor, &Hp[j+np], &rp[j]);
            mclBnG1_add(&R_vec[i], &R_vec[i], &g1Factor);
        }

        mclBnFr_setByCSPRNG(&xp_vec[i]);
        mclBnFr_inv(&frFactor, &xp_vec[i]);

        for (int j = 0; j < np; j++)
        {
            mclBnG1_mul(&g1Factor, &Gp[j], &frFactor);
            mclBnG1_mul(&Gp[j], &Gp[j+np], &xp_vec[i]);
            mclBnG1_add(&Gp[j], &Gp[j], &g1Factor);

            mclBnG1_mul(&Hp[j], &Hp[j], &xp_vec[i]);
            mclBnG1_mul(&g1Factor, &Hp[j+np], &frFactor);
            mclBnG1_add(&Hp[j], &Hp[j], &g1Factor);

            mclBnFr_mul(&lp[j], &lp[j], &xp_vec[i]);
            mclBnFr_mul(&frFactor2, &lp[j+np], &frFactor);
            mclBnFr_add(&lp[j], &lp[j], &frFactor2);

            mclBnFr_mul(&frFactor2, &rp[j], &frFactor);
            mclBnFr_mul(&rp[j], &rp[j+np], &xp_vec[i]);
            mclBnFr_add(&rp[j], &rp[j], &frFactor2);
        }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    printf("\033[1;32m[SUCCESS] :\033[0m Bulletproof created in ");
    printf("%fs\n", elapsed);

    // VERIFIER
    clock_gettime(CLOCK_MONOTONIC, &begin);

    mclBnFr delta_yz;
    mclBnFr_clear(&frFactor);
    mclBnFr_clear(&frFactor2);

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_add(&frFactor, &frFactor, &y_vec[i]);
    }

    for (int i = 0; i < Nb; i++)
    {
        mclBnFr_add(&frFactor2, &frFactor2, &two_vec[i]);
    }

    mclBnFr_sub(&delta_yz, &z, &z2);
    mclBnFr_mul(&delta_yz, &delta_yz, &frFactor);

    mclBnFr_mul(&frFactor, &one, &z2);
    for (int i = 1; i < Mc + 1; i++)
    {
        mclBnFr_mul(&frFactor, &frFactor, &z);
        mclBnFr_mul(&frFactor3, &frFactor, &frFactor2);
        mclBnFr_add(&frFactor4, &frFactor4, &frFactor3);
    }
    
    mclBnFr_sub(&delta_yz, &delta_yz, &frFactor4);

    mclBnFr z_vec[Mc];
    mclBnFr_setInt(&z_vec[0], 1);
    for (int i = 1; i < Mc; i++)
    {
        mclBnFr_mul(&z_vec[i], &z_vec[i-1], &z);
    }

    mclBnG1 CL, CR;
    mclBnG1_mul(&g1Factor, &Gb, &t_inner);
    mclBnG1_mul(&CL, &Hb, &tx);
    mclBnG1_add(&CL, &CL, &g1Factor);
    mclBnG1_clear(&CR);

    for (int i = 0; i < Mc; i++)
    {
        mclBnFr_mul(&frFactor, &z_vec[i], &z2);
        mclBnG1_mul(&g1Factor, &V[i], &frFactor);
        mclBnG1_add(&CR, &CR, &g1Factor);
    }
    
    mclBnG1_mul(&g1Factor, &Gb, &delta_yz);
    mclBnG1_add(&CR, &CR, &g1Factor);
    mclBnG1_mul(&g1Factor, &T1, &x);
    mclBnG1_add(&CR, &CR, &g1Factor);
    mclBnG1_mul(&g1Factor, &T2, &x);
    mclBnG1_mul(&g1Factor, &g1Factor, &x);
    mclBnG1_add(&CR, &CR, &g1Factor);   

    mclBnG1_getStr(buff1, sizeof(buff1), &CL, 10);
    mclBnG1_getStr(buff2, sizeof(buff2), &CR, 10);

    int cond1 = 0;
    int cond2 = 0;

    if (!strcmp(buff1, buff2)) cond1 = 1;

    mclBnFr s_vec[Nb*Mc];

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_setInt(&s_vec[i], 1);
        for (int j = 0; j < logN; j++)
        {
            int bin = i >> (logN - 1 - j);
            if (bin & 1) mclBnFr_mul(&s_vec[i], &s_vec[i], &xp_vec[j]);
            else
            {
                mclBnFr_inv(&frFactor, &xp_vec[j]);
                mclBnFr_mul(&s_vec[i], &s_vec[i], &frFactor);
            }
        }
    }

    mclBnG1 LHS, RHS;

    mclBnFr_mul(&frFactor, &lp[0], &rp[0]);
    mclBnG1_mul(&LHS, &Ub, &frFactor);

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&frFactor, &s_vec[i], &lp[0]);
        mclBnG1_mul(&g1Factor, &G[i], &frFactor);
        mclBnG1_add(&LHS, &LHS, &g1Factor);

        mclBnFr_inv(&frFactor, &s_vec[i]);
        mclBnFr_mul(&frFactor, &frFactor, &rp[0]);
        mclBnG1_mul(&g1Factor, &Hs[i], &frFactor);
        mclBnG1_add(&LHS, &LHS, &g1Factor);
    }

    mclBnG1_clear(&RHS);
    mclBnG1_add(&RHS, &RHS, &Pp);

    for (int i = 0; i < logN; i++)
    {
        mclBnFr_mul(&frFactor, &xp_vec[i], &xp_vec[i]);
        mclBnG1_mul(&g1Factor, &L_vec[i], &frFactor);
        mclBnG1_add(&RHS, &RHS, &g1Factor);

        mclBnFr_mul(&frFactor, &xp_vec[i], &xp_vec[i]);
        mclBnFr_inv(&frFactor, &frFactor);
        mclBnG1_mul(&g1Factor, &R_vec[i], &frFactor);
        mclBnG1_add(&RHS, &RHS, &g1Factor);
    }

    mclBnG1_getStr(buff1, sizeof(buff1), &LHS, 10);
    mclBnG1_getStr(buff2, sizeof(buff2), &RHS, 10);

    if (!strcmp(buff1, buff2)) cond2 = 1;

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    if ((cond1) && (cond2))
    {
        printf("\033[1;32m[SUCCESS] :\033[0m Bulletproof verified in ");
        printf("%fs\n", elapsed);
    } 
    else printf("\033[1;31m[FAIL] :\033[0m Bulletproof INCORRECT.");
}