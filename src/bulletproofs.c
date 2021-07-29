mclBnG1 Gb, Hb, Ub;
mclBnG1 *G, *H;
mclBnG1 *V;
mclBnG1 A, S;
mclBnFr y, z;
mclBnFr x, t_inner, tx, mu;
mclBnG1 T1, T2;
mclBnG1 *Hp;
mclBnFr xp;
mclBnFr *xp_vec;
mclBnG1 *L_vec, *R_vec;
mclBnFr *l, *r;
mclBnFr *two_vec;
mclBnFr *gammas;

int userGammas = 0;

void bulletproof_prove(unsigned char *si[])
{
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    mclBnG1 Gen;
    mclBnG1_setStr(&Gen, GGEN, strlen(GGEN), 10);

    mclBnFr rnd;
    mclBnG1 Ubn;

    if (!userGammas)
    {
        for (int i = 0; i < Mc; i++)
        {
            mclBnFr_setByCSPRNG(&gammas[i]);
        }
    }
    
    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Gb, &Gen, &rnd);

    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Hb, &Gen, &rnd);

    mclBnFr_setByCSPRNG(&rnd);
    mclBnG1_mul(&Ub, &Gen, &rnd);

    mclBnFr frFactor, frFactor2, frFactor3, frFactor4;

    mclBnG1_clear(&G[0]);
    mclBnG1_clear(&H[0]); 

    for (int i = 1; i < Nb*Mc; i++)
    {
        mclBnG1_add(&G[i], &G[i-1], &Gb);
        mclBnG1_add(&H[i], &H[i-1], &Hb);
    }

    mclBnFr aL[Nb*Mc], aR[Nb*Mc];

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_clear(&aL[i]);
    }

    mclBnFr v[Mc];
    mclBnG1 g1Factor;
    char buff[2048];

    for (int j = 0; j < Mc; j++)
    {
        mclBnFr_setStr(&v[j], si[j], strlen(si[j]), 10);
        mclBnFr_getStr(buff, sizeof(buff), &v[j], 2);

        for (int i = 0; i < strlen(buff); i++)
        {
            if(buff[strlen(buff) - 1 - i] == '1') mclBnFr_setInt(&aL[i + (j*Nb)], 1);
        }
    }

    mclBnFr one;
    mclBnFr_setInt(&one, 1);
    mclBnFr alpha;
    mclBnFr_setByCSPRNG(&alpha);

    mclBnFr sL[Nb*Mc], sR[Nb*Mc];

    mclBnFr rho;
    mclBnFr_setByCSPRNG(&rho);

    mclBnG1_mul(&A, &Hb, &alpha);
    mclBnG1_mul(&S, &Hb, &rho);

    #pragma omp parallel for
    for (int i = 0; i < Mc; i++)
    {
        // compute commitment V
        mclBnG1 g1Tmp;
        mclBnG1_mul(&g1Tmp, &Gb, &v[i]);
        mclBnG1_mul(&V[i], &Hb, &gammas[i]);
        mclBnG1_add(&V[i], &V[i], &g1Tmp);
    }

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_sub(&aR[i], &aL[i], &one);
        
        mclBnFr_setByCSPRNG(&sL[i]);
        mclBnFr_setByCSPRNG(&sR[i]);
    }

    // compute commitment A
    for (int i = 0; i < Nb*Mc; i++)
    {
        if (mclBnFr_isOne(&aL[i])) mclBnG1_add(&A, &A, &G[i]);
    }

    mclBnG1 A_chunk;
    mult_exp(&A_chunk, H, aR, Nb*Mc);
    mclBnG1_add(&A, &A, &A_chunk);

    // compute commitment S
    mclBnG1 S_chunk;
    mult_exp(&S_chunk, G, sL, Nb*Mc);
    mclBnG1_add(&S, &S, &S_chunk);
    mult_exp(&S_chunk, H, sR, Nb*Mc);
    mclBnG1_add(&S, &S, &S_chunk);

    transcript_add_G1(&A);
    transcript_add_G1(&S);

    mclBnFr z2;
    transcript_hash(&y);
    transcript_add_Fr(&y);
    transcript_hash(&z);
    transcript_add_Fr(&z);

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
    mclBnFr_setByCSPRNG(&tau1);
    mclBnFr_setByCSPRNG(&tau2);

    mclBnG1_mul(&g1Factor, &Gb, &t1);
    mclBnG1_mul(&T1, &Hb, &tau1);
    mclBnG1_add(&T1, &T1, &g1Factor);

    mclBnG1_mul(&g1Factor, &Gb, &t2);
    mclBnG1_mul(&T2, &Hb, &tau2);
    mclBnG1_add(&T2, &T2, &g1Factor);

    transcript_add_G1(&T1);
    transcript_add_G1(&T2);

    transcript_hash(&x);
    transcript_add_Fr(&x);

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
        mclBnFr_mul(&frFactor, &frFactor2, &gammas[i-1]);
        mclBnFr_add(&frFactor3, &frFactor3, &frFactor);
    }

    mclBnFr_add(&tx, &tx, &frFactor3);

    mclBnFr_mul(&mu, &rho, &x);
    mclBnFr_add(&mu, &mu, &alpha);

    // inner product
    transcript_hash(&xp);
    transcript_add_Fr(&xp);
    
    mclBnG1_mul(&Ubn, &Ub, &xp);

    mclBnG1 Gp[Nb*Mc];

    #pragma omp parallel for
    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnG1_mul(&Gp[i], &G[i], &one);

        mclBnFr frTmp;
        mclBnFr_inv(&frTmp, &y_vec[i]);
        mclBnG1_mul(&Hp[i], &H[i], &frTmp);
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
            mclBnFr_mul(&frFactor, &l[j], &r[j+np]);
            mclBnFr_add(&cl, &cl, &frFactor);
            mclBnFr_mul(&frFactor, &l[j+np], &r[j]);
            mclBnFr_add(&cr, &cr, &frFactor);
        }

        mclBnG1 L, R;
        mclBnG1_mul(&L, &Ubn, &cl);
        mclBnG1_mul(&R, &Ubn, &cr);

        mclBnG1 L_chunk;
        mult_exp(&L_chunk, Gp+np, l, np);
        mclBnG1_add(&L, &L, &L_chunk);
        mult_exp(&L_chunk, Hp, r+np, np);
        mclBnG1_add(&L, &L, &L_chunk);

        mclBnG1 R_chunk;
        mult_exp(&R_chunk, Gp, l+np, np);
        mclBnG1_add(&R, &R, &R_chunk);
        mult_exp(&R_chunk, Hp+np, r, np);
        mclBnG1_add(&R, &R, &R_chunk);

        transcript_add_G1(&L);
        transcript_add_G1(&R);

        transcript_hash(&xp);
        transcript_add_Fr(&xp);

        mclBnFr_inv(&frFactor, &xp);

        #pragma omp parallel for
        for (int j = 0; j < np; j++)
        {
            mclBnG1 g1Tmp;
            mclBnFr frTmp;
            mclBnG1_mul(&g1Tmp, &Gp[j], &frFactor);
            mclBnG1_mul(&Gp[j], &Gp[j+np], &xp);
            mclBnG1_add(&Gp[j], &Gp[j], &g1Tmp);

            mclBnG1_mul(&Hp[j], &Hp[j], &xp);
            mclBnG1_mul(&g1Tmp, &Hp[j+np], &frFactor);
            mclBnG1_add(&Hp[j], &Hp[j], &g1Tmp);

            mclBnFr_mul(&l[j], &l[j], &xp);
            mclBnFr_mul(&frTmp, &l[j+np], &frFactor);
            mclBnFr_add(&l[j], &l[j], &frTmp);

            mclBnFr_mul(&frTmp, &r[j], &frFactor);
            mclBnFr_mul(&r[j], &r[j+np], &xp);
            mclBnFr_add(&r[j], &r[j], &frTmp);
        }
    }
    
    bulletproof_save();

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    if (bench)
    {
        printf("\033[1;32m[SUCCESS] :\033[0m Bulletproof created in ");
        printf("%fs\n", elapsed);
    }
}

int bulletproof_verify()
{
    struct timespec begin, end;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &begin);

    mclBnFr frFactor, frFactor2, frFactor3, frFactor4;

    mclBnFr one;
    mclBnFr_setInt(&one, 1);
    mclBnG1 g1Factor;

    mclBnFr z2;
    mclBnFr y_vec[Nb*Mc];
    mclBnG1 P, Pp;

    bulletproof_read();

    mclBnG1_clear(&G[0]);
    mclBnG1_clear(&H[0]); 

    for (int i = 1; i < Nb*Mc; i++)
    {
        mclBnG1_add(&G[i], &G[i-1], &Gb);
        mclBnG1_add(&H[i], &H[i-1], &Hb);
    }

    mclBnFr_mul(&z2, &z, &z);
    mclBnFr_setInt(&y_vec[0], 1);
    mclBnFr_mul(&y_vec[1], &y_vec[0], &y);

    for (int i = 2; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&y_vec[i], &y_vec[i-1], &y);
    }

    mclBnG1_mul(&P, &S, &x);
    mclBnG1_add(&P, &P, &A);

    mclBnFr_clear(&frFactor3);
    mclBnFr_mul(&frFactor3, &one, &z);

    #pragma omp parallel for
    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr frTmp;
        mclBnFr_inv(&frTmp, &y_vec[i]);
        mclBnG1_mul(&Hp[i], &H[i], &frTmp);
    }

    mclBnFr frFactor_vec[Nb*Mc];
    mclBnFr frFactor_vec2[Nb*Mc];
    mclBnG1 P_chunk;

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&frFactor, &z, &y_vec[i]);
        if (i % Nb == 0) mclBnFr_mul(&frFactor3, &frFactor3, &z);
        mclBnFr_mul(&frFactor2, &frFactor3, &two_vec[i % Nb]);
        mclBnFr_add(&frFactor_vec[i], &frFactor, &frFactor2);

        mclBnFr_neg(&frFactor_vec2[i], &z);
    }

    mult_exp(&P_chunk, Hp, frFactor_vec, Nb*Mc);
    mclBnG1_add(&P, &P, &P_chunk);
    mult_exp(&P_chunk, G, frFactor_vec2, Nb*Mc);
    mclBnG1_add(&P, &P, &P_chunk);

    mclBnFr_neg(&frFactor, &mu);
    mclBnG1_mul(&g1Factor, &Hb, &frFactor);
    mclBnG1_add(&Pp, &P, &g1Factor);

    mclBnFr_mul(&frFactor, &xp, &t_inner);
    mclBnG1_mul(&g1Factor, &Ub, &frFactor);
    mclBnG1_add(&Pp, &Pp, &g1Factor);

    mclBnG1_mul(&Ub, &Ub, &xp);

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
    mclBnFr_clear(&frFactor4);

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
        mclBnFr_mul(&frFactor_vec[i], &z_vec[i], &z2);
    }

    mult_exp(&CR, V, frFactor_vec, Mc);
    
    mclBnG1_mul(&g1Factor, &Gb, &delta_yz);
    mclBnG1_add(&CR, &CR, &g1Factor);
    mclBnG1_mul(&g1Factor, &T1, &x);
    mclBnG1_add(&CR, &CR, &g1Factor);
    mclBnG1_mul(&g1Factor, &T2, &x);
    mclBnG1_mul(&g1Factor, &g1Factor, &x);
    mclBnG1_add(&CR, &CR, &g1Factor);   

    int cond1 = mclBnG1_isEqual(&CL, &CR);

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
    mclBnG1 LHS_chunk, RHS_chunk;

    mclBnFr_mul(&frFactor, &l[0], &r[0]);
    mclBnG1_mul(&LHS, &Ub, &frFactor);

    for (int i = 0; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&frFactor_vec[i], &s_vec[i], &l[0]);
        mclBnFr_inv(&frFactor, &s_vec[i]);
        mclBnFr_mul(&frFactor_vec2[i], &frFactor, &r[0]);
    }

    mult_exp(&LHS_chunk, G, frFactor_vec, Nb*Mc);
    mclBnG1_add(&LHS, &LHS, &LHS_chunk);
    mult_exp(&LHS_chunk, Hp, frFactor_vec2, Nb*Mc);
    mclBnG1_add(&LHS, &LHS, &LHS_chunk);

    mclBnG1_clear(&RHS);
    mclBnG1_add(&RHS, &RHS, &Pp);

    for (int i = 0; i < logN; i++)
    {
        mclBnFr_mul(&frFactor_vec[i], &xp_vec[i], &xp_vec[i]);
        mclBnFr_inv(&frFactor_vec2[i], &frFactor_vec[i]);
    }

    mult_exp(&RHS_chunk, L_vec, frFactor_vec, logN);
    mclBnG1_add(&RHS, &RHS, &RHS_chunk);
    mult_exp(&RHS_chunk, R_vec, frFactor_vec2, logN);
    mclBnG1_add(&RHS, &RHS, &RHS_chunk);

    int cond2 = mclBnG1_isEqual(&LHS, &RHS);

    clock_gettime(CLOCK_MONOTONIC, &end);
    elapsed = (end.tv_sec - begin.tv_sec);
    elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;

    if ((cond1) && (cond2))
    {
        if (bench)
        {
            printf("\033[1;32m[SUCCESS] :\033[0m Bulletproof verified in ");
            printf("%fs\n", elapsed);
        }

        return 1;
    } 
    else
    {
        if (bench) printf("\033[1;31m[FAIL] :\033[0m Bulletproof INCORRECT."); 
        return 0;
    }
}

void bulletproof_save()
{
    char buff[2048];

    FILE *fbp;
    fbp = fopen("data/bulletproof.params", "w");

    mclBnG1_getStr(buff, sizeof(buff), &Gb, 10);
    fprintf(fbp, "%s\n", buff);
    mclBnG1_getStr(buff, sizeof(buff), &Hb, 10);
    fprintf(fbp, "%s\n", buff);
    mclBnG1_getStr(buff, sizeof(buff), &Ub, 10);
    fprintf(fbp, "%s\n", buff);

    for (int i = 0; i < Mc; i++)
    {
        mclBnG1_getStr(buff, sizeof(buff), &V[i], 10);
        fprintf(fbp, "%s\n", buff);
    }

    mclBnFr_getStr(buff, sizeof(buff), &tx, 10);
    fprintf(fbp, "%s\n", buff);
    mclBnFr_getStr(buff, sizeof(buff), &t_inner, 10);
    fprintf(fbp, "%s\n", buff);
    mclBnFr_getStr(buff, sizeof(buff), &mu, 10);
    fprintf(fbp, "%s\n", buff);

    mclBnFr_getStr(buff, sizeof(buff), &l[0], 10);
    fprintf(fbp, "%s\n", buff);
    mclBnFr_getStr(buff, sizeof(buff), &r[0], 10);
    fprintf(fbp, "%s\n", buff);

    fprintf(fbp, "%s\n", transcript);

    fclose(fbp);
}

void bulletproof_read()
{
    char buff[2048];

    FILE *fbp;
    fbp = fopen("data/bulletproof.params", "r");

    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&Gb, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&Hb, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&Ub, buff, strlen(buff), 10);

    for (int i = 0; i < Mc; i++)
    {
        fgets(buff, sizeof buff, fbp);
        mclBnG1_setStr(&V[i], buff, strlen(buff), 10);
    }

    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&tx, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&t_inner, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&mu, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&l[0], buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&r[0], buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&A, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&S, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&y, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&z, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&T1, buff, strlen(buff), 10);
    fgets(buff, sizeof buff, fbp);
    mclBnG1_setStr(&T2, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&x, buff, strlen(buff), 10);

    fgets(buff, sizeof buff, fbp);
    mclBnFr_setStr(&xp, buff, strlen(buff), 10);

    for (int i = 0; i < logN; i++)
    {
        fgets(buff, sizeof buff, fbp);
        mclBnG1_setStr(&L_vec[i], buff, strlen(buff), 10);
        fgets(buff, sizeof buff, fbp);
        mclBnG1_setStr(&R_vec[i], buff, strlen(buff), 10);  
        fgets(buff, sizeof buff, fbp);
        mclBnFr_setStr(&xp_vec[i], buff, strlen(buff), 10);      
    }

    fclose(fbp);
}

static inline void bulletproof_init(int Nb_set, int Mc_set)
{
    Nb = Nb_set;
    Mc = Mc_set;

    struct stat st = {0};
    if (stat("data", &st) == -1) mkdir("data", 0700);

    mclBn_init(USEDCURVE, MCLBN_COMPILED_TIME_VAR);

    float log_up = log(Nb*Mc);
    float log_down = log(2);
    logN = log_up / log_down;

    transcript = (char *) malloc(1024 * logN * sizeof(char));

    G = (mclBnG1*) malloc((Nb*Mc) * sizeof(mclBnG1));
    H = (mclBnG1*) malloc((Nb*Mc) * sizeof(mclBnG1));
    V = (mclBnG1*) malloc((Mc) * sizeof(mclBnG1));
    Hp = (mclBnG1*) malloc((Nb*Mc) * sizeof(mclBnG1));
    xp_vec = (mclBnFr*) malloc((logN) * sizeof(mclBnFr));
    L_vec = (mclBnG1*) malloc((logN) * sizeof(mclBnG1));
    R_vec = (mclBnG1*) malloc((logN) * sizeof(mclBnG1));
    l = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    r = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    two_vec = (mclBnFr*) malloc((Nb*Mc) * sizeof(mclBnFr));
    gammas = (mclBnFr*) malloc((Mc) * sizeof(mclBnFr));

    mclBnFr_setInt(&two_vec[0], 1);
    mclBnFr_setInt(&two_vec[1], 2);

    for (int i = 2; i < Nb*Mc; i++)
    {
        mclBnFr_mul(&two_vec[i], &two_vec[i-1], &two_vec[1]);
    }
}

static inline void bulletproof_get_context(context *ctx)
{
    ctx->V = V;
    ctx->G = Gb;
    ctx->H = Hb;
}

static inline void bulletproof_get_gammas(context *ctx)
{
    ctx->gammas = gammas;
}

static inline void bulletproof_user_gammas(int val)
{
    userGammas = val;
}