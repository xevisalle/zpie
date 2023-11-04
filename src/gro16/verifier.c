
int verify(mclBnG1 *piA, mclBnG2 *piB2, mclBnG1 *piC, mpz_t u[(nPublic + nConst)], verifying_key vk)
{
    mclBnG1 factorG1;
    mclBnFr frFactor;
    mclBnGT pairing1, pairing2, pairing3, factorGT;
    mclBnG1 Vu;

    mclBnG1_clear(&Vu);
    for (int i = (nPublic); i--;)
    {
        // Vu = Vu + u[i] * s1.vk[i]
        mpz_to_fr(&frFactor, &u[i]);
        mclBnG1_mul(&factorG1, &vk.vk1[i+nConst], &frFactor);
        mclBnG1_add(&Vu, &Vu, &factorG1);
    }

    for (int i = (nConst); i--;)
    {
        // Vu = Vu + u[i] * s1.vk[i]
        mpz_to_fr(&frFactor, &vk.constants[i]);
        mclBnG1_mul(&factorG1, &vk.vk1[i], &frFactor);
        mclBnG1_add(&Vu, &Vu, &factorG1);
    }

    log_message("Computing e(piA, piB2), e(Vu, gamma), e(piC, delta)...");

    #pragma omp parallel num_threads(3)
    {
        switch (get_thread())
        {
            case 0: mclBn_pairing(&pairing1, piA, piB2); break;
            case 1: mclBn_pairing(&pairing2, &Vu, &vk.gamma2); break;
            case 2: mclBn_pairing(&pairing3, piC, &vk.delta2); break;
            case 99:
                mclBn_pairing(&pairing1, piA, piB2);
                mclBn_pairing(&pairing2, &Vu, &vk.gamma2);
                mclBn_pairing(&pairing3, piC, &vk.delta2);
                break;
        }
    }

    log_state(1);
    char buff[2048];
    if (logs)
    {
        mclBnGT_getStr(buff, sizeof(buff), &pairing1, 10);
        printf("e(piA, piB2) = %s\n", buff);
        mclBnGT_getStr(buff, sizeof(buff), &pairing2, 10);
        printf("e(Vu, gamma) = %s\n", buff);
        mclBnGT_getStr(buff, sizeof(buff), &pairing3, 10);
        printf("e(piC, delta) = %s\n", buff);
    }

    log_message("Computing e(alpha, beta) * e(Vu, gamma) * e(piC, delta)...");
    mclBnGT_mul(&factorGT, &vk.alphabetaT, &pairing2);
    mclBnGT_mul(&factorGT, &factorGT, &pairing3);
    log_state(1);
    if (logs)
    {
        mclBnGT_getStr(buff, sizeof(buff), &factorGT, 10);
        printf("e(alpha, beta) * e(Vu, gamma) * e(piC, delta) = %s\n", buff);
    }

    log_message("e(piA, piB2) = e(alpha, beta) * e(Vu, gamma) * e(piC, delta)???");
    
    int verified = mclBnGT_isEqual(&pairing1, &factorGT);

    if (mclBnGT_isOne(&pairing2)) verified = 0;
    
    if (verified) log_state(1);
    else log_state(0);
    
    return verified;
}