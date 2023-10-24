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

    CU_basic_run_tests();

    CU_cleanup_registry();
    return CU_get_error();
}