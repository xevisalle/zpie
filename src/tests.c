#include "zpie.h"

#include "../circuits/mimc.c"
#include "../circuits/eddsa.c"

void test_single_constraint()
{
    element out;
    init_public(&out);

    element a, b;
    init(&a);
    init(&b);

    input(&a, "1234");
    input(&b, "5678");

    mul(&out, &a, &b);
}

void test_eddsa_verification()
{
    point B;
    B.x = "5299619240641551281634865583518297030282874472190772894086521144482721001553";
    B.y = "16950150798460657717958625567821834550301663161624707787222815936182638968203";
    
    eddsa_signature edsig;
    edsig.R.x = "1262948111445225057373438194818763405700457487429548371463214326190311895864"; 
    edsig.R.y = "12533500305127747239777484416561675628195562065959201739446841668623540883587";
    edsig.S = "2674591880888862378688383832785447197125897205360861957116147165712709455207";

    point A;
    A.x = "21629779320182474195265732521833299809982444552305142529409236301104997786342";
    A.y = "9011812445381030664142622066218331845140881847034934166630871421746105699091";

    char *msg = "1234";

    verify_eddsa(edsig, B, A, msg);
}

void test_mimc_hash()
{
    element h, x_in, k;

    init_public(&h);
    init(&x_in);
    init(&k);
    
    input(&x_in, "1234");
    input(&k, "112233445566");

    mimc7(&h, &x_in, &k);
}

void test_setup(void)
{
    test_no_rand = 1;
    setup_keys keys = perform_setup(&test_single_constraint); 

    char* pk_bytes = serialize_pk(&keys.pk);
    char* vk_bytes = serialize_vk(&keys.vk);

    BYTE hash_bytes[SHA256_BLOCK_SIZE];
    SHA256_CTX ctx;

    sha256_init(&ctx);
    sha256_update(&ctx, pk_bytes, strlen(pk_bytes));
    sha256_final(&ctx, hash_bytes);

    CU_ASSERT(!strcmp(to_hex(hash_bytes, sizeof hash_bytes), "b8a812b4c6576d343f0c269157a0020ca63e808244f2e0f75b9940797502d4fa"));

    sha256_init(&ctx);
    sha256_update(&ctx, vk_bytes, strlen(vk_bytes));
    sha256_final(&ctx, hash_bytes);

    CU_ASSERT(!strcmp(to_hex(hash_bytes, sizeof hash_bytes), "dac6dfe723f1874422a4235e38a6f4ac4bb1b18716c90babe3d4a97f189ac15e"));

    test_no_rand = 0;
}

void test_prover(void)
{
    test_no_rand = 1;
    setup_keys keys = perform_setup(&test_single_constraint); 
    proof p = generate_proof(&test_single_constraint, keys.pk);

    const char *piAstr = "1 13398732126763033363928255770670403609664455533535809960659793057603927642327 14567332642717250669329472598965177550050834309459245026995104363234319745805";
    const char *piB2str = "1 9513526328373247288214002967710658327692956864193416721895179753121227228903 17320346092699268035923233491595138958007151833266586455159840335219170425243 8079768110185479532548096263199181437927983909022782182442306192699700743609 19381997603489315175356927627025590277145986935796790438444340629346184509934";
    const char *piCstr = "1 6751941069502688487334371509286578067074020223942252322110100779175835131489 10460091663676025417104943726359531715081933829156881875323036992094404259688";

    mclBnG1 piA, piC;
    mclBnG2 piB2;

    mclBnG1_setStr(&piA, piAstr, strlen(piAstr), 10);
    mclBnG2_setStr(&piB2, piB2str, strlen(piB2str), 10);
    mclBnG1_setStr(&piC, piCstr, strlen(piCstr), 10);

    CU_ASSERT(mclBnG1_isEqual(&piA, &p.piA));
    CU_ASSERT(mclBnG2_isEqual(&piB2, &p.piB2));
    CU_ASSERT(mclBnG1_isEqual(&piC, &p.piC));

    test_no_rand = 0;
}

void test_full_circuits(void)
{
    setup_keys keys_sc = perform_setup(&test_single_constraint); 
    setup_keys keys_mh = perform_setup(&test_mimc_hash);  
    setup_keys keys_ev = perform_setup(&test_eddsa_verification); 

    proof p_sc = generate_proof(&test_single_constraint, keys_sc.pk);
    proof p_mh = generate_proof(&test_mimc_hash, keys_mh.pk);
    proof p_ev = generate_proof(&test_eddsa_verification, keys_ev.pk);
    
    CU_ASSERT(verify_proof(&test_single_constraint, p_sc, keys_sc.vk));
    CU_ASSERT(verify_proof(&test_mimc_hash, p_mh, keys_mh.vk));
    CU_ASSERT(verify_proof(&test_eddsa_verification, p_ev, keys_ev.vk));
}

//TODO: fix this
void test_bulletproofs(void)
{
    // we init the bulletproofs module, for 2 aggregated proofs of 64 bits
    //bulletproof_init(64, 2);

    // we set some values to prove knowledge of, and compute the proof (../data/bulletproof.params)
    //unsigned char *si[] = {"1234", "5678"};
    //bulletproof_prove(si);

    // we verify the bulletproof (../data/bulletproof.params)
    //if(bulletproof_verify()) printf("Bulletproof verified.\n");
    //else printf("Bulletproof cannot be verified.\n");
}

int main()
{
    CU_pSuite suite = NULL;
    if (CUE_SUCCESS != CU_initialize_registry()) return CU_get_error();
    suite = CU_add_suite("Test Suite", init_suite, clean_suite);

    if ((NULL == suite) || (NULL == CU_add_test(suite, "\n\nFull Circuits Testing\n\n", test_full_circuits)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    if ((NULL == suite) || (NULL == CU_add_test(suite, "\n\nConstraint System Testing\n\n", test_constraint_system)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    if ((NULL == suite) || (NULL == CU_add_test(suite, "\n\nProver Testing\n\n", test_prover)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    if ((NULL == suite) || (NULL == CU_add_test(suite, "\n\nSetup Testing\n\n", test_setup)))
    {
        CU_cleanup_registry();
        return CU_get_error();
    }

    CU_basic_run_tests();

    CU_cleanup_registry();
    return CU_get_error();
}