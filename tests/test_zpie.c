#include <criterion/criterion.h>
#include <zpie.h>
#include <zpie_internal.h>

#include "../gadgets/eddsa.c"
#include "../gadgets/mimc.c"

void test_single_constraint()
{
    zpie_element out;
    zpie_init_public(&out);

    zpie_element a, b;
    zpie_init(&a);
    zpie_init(&b);

    zpie_input(&a, "1234");
    zpie_input(&b, "5678");

    zpie_mul(&out, &a, &b);
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

    char* msg = "1234";

    verify_eddsa(edsig, B, A, msg);
}

void test_mimc_hash()
{
    zpie_element h, x_in, k;

    zpie_init_public(&h);
    zpie_init(&x_in);
    zpie_init(&k);

    zpie_input(&x_in, "1234");
    zpie_input(&k, "112233445566");

    mimc7(&h, &x_in, &k);
}

Test(zpie, prover)
{
    test_no_rand = 1;
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_single_constraint);
    zpie_proof p;
    zpie_generate_proof(&p, &test_single_constraint, &keys.pk);

    const char* piAstr =
        "1 13398732126763033363928255770670403609664455533535809960659793057603927642327 "
        "14567332642717250669329472598965177550050834309459245026995104363234319745805";
    const char* piB2str =
        "1 9513526328373247288214002967710658327692956864193416721895179753121227228903 "
        "17320346092699268035923233491595138958007151833266586455159840335219170425243 "
        "8079768110185479532548096263199181437927983909022782182442306192699700743609 "
        "19381997603489315175356927627025590277145986935796790438444340629346184509934";
    const char* piCstr =
        "1 6751941069502688487334371509286578067074020223942252322110100779175835131489 "
        "10460091663676025417104943726359531715081933829156881875323036992094404259688";

    mclBnG1 piA, piC;
    mclBnG2 piB2;

    mclBnG1_setStr(&piA, piAstr, strlen(piAstr), 10);
    mclBnG2_setStr(&piB2, piB2str, strlen(piB2str), 10);
    mclBnG1_setStr(&piC, piCstr, strlen(piCstr), 10);

    cr_assert(mclBnG1_isEqual(&piA, &p.piA));
    cr_assert(mclBnG2_isEqual(&piB2, &p.piB2));
    cr_assert(mclBnG1_isEqual(&piC, &p.piC));

    test_no_rand = 0;
}

Test(zpie, full_circuits)
{
    zpie_setup_keys keys_sc;
    zpie_perform_setup(&keys_sc, &test_single_constraint);
    zpie_setup_keys keys_mh;
    zpie_perform_setup(&keys_mh, &test_mimc_hash);
    zpie_setup_keys keys_ev;
    zpie_perform_setup(&keys_ev, &test_eddsa_verification);

    zpie_proof p_sc;
    zpie_generate_proof(&p_sc, &test_single_constraint, &keys_sc.pk);
    zpie_proof p_mh;
    zpie_generate_proof(&p_mh, &test_mimc_hash, &keys_mh.pk);
    zpie_proof p_ev;
    zpie_generate_proof(&p_ev, &test_eddsa_verification, &keys_ev.pk);

    cr_assert(zpie_verify_proof(&test_single_constraint, &p_sc, &keys_sc.vk));
    cr_assert(zpie_verify_proof(&test_mimc_hash, &p_mh, &keys_mh.vk));
    cr_assert(zpie_verify_proof(&test_eddsa_verification, &p_ev, &keys_ev.vk));
}

void test_full_api()
{
    zpie_element e_mul, e_addmul, e_add3mul, e_addmuladd;
    zpie_init(&e_mul);
    zpie_init(&e_addmul);
    zpie_init(&e_add3mul);
    zpie_init(&e_addmuladd);

    zpie_element a, b;
    zpie_init(&a);
    zpie_init(&b);

    zpie_input(&a, "5");
    zpie_input(&b, "10");

    zpie_mul(&e_mul, &a, &b);
    zpie_addmul(&e_addmul, &a, &b, &b);
    zpie_add3mul(&e_add3mul, &a, &a, &a, &b);
    zpie_addmuladd(&e_addmuladd, &a, &a, &b, &b);
}

void test_all_constraints()
{
    zpie_element a, b, c;
    zpie_init_public(&a);
    zpie_init(&b);
    zpie_init(&c);

    zpie_input(&a, "3");
    zpie_input(&b, "7");
    zpie_input(&c, "5");

    // zpie_mul: out = a * b = 21
    zpie_element e_mul;
    zpie_init(&e_mul);
    zpie_mul(&e_mul, &a, &b);

    // zpie_addmul: out = (a + b) * c = 50
    zpie_element e_addmul;
    zpie_init(&e_addmul);
    zpie_addmul(&e_addmul, &a, &b, &c);

    // zpie_add3mul: out = (a + b + c) * a = 45
    zpie_element e_add3mul;
    zpie_init(&e_add3mul);
    zpie_add3mul(&e_add3mul, &a, &b, &c, &a);

    // zpie_addmuladd: out = (a + b) * (c + a) = 80
    zpie_element e_addmuladd;
    zpie_init(&e_addmuladd);
    zpie_addmuladd(&e_addmuladd, &a, &b, &c, &a);

    // zpie_add3muladd3: out = (a + b + c) * (a + b + c) = 225
    zpie_element e_add3muladd3;
    zpie_init(&e_add3muladd3);
    zpie_add3muladd3(&e_add3muladd3, &a, &b, &c, &a, &b, &c);

    // zpie_mul_constants: out = (2*a) * (3*b) = 126
    zpie_element e_mulc;
    zpie_init(&e_mulc);
    int lc = 2, rc = 3;
    zpie_mul_constants(&e_mulc, &lc, &a, &rc, &b);

    // zpie_addmul_constants: out = (2*a + 3*b) * (4*c) = 540
    zpie_element e_addmulc;
    zpie_init(&e_addmulc);
    int lc1 = 2, lc2 = 3, rc2 = 4;
    zpie_addmul_constants(&e_addmulc, &lc1, &a, &lc2, &b, &rc2, &c);

    // zpie_mul_big_constants: out = (10*a) * (20*b) = 4200
    zpie_element e_bigc;
    zpie_init(&e_bigc);
    mclBnFr big_lc, big_rc;
    mclBnFr_setInt(&big_lc, 10);
    mclBnFr_setInt(&big_rc, 20);
    zpie_mul_big_constants(&e_bigc, &big_lc, &a, &big_rc, &b);

    // zpie_addsmul: out = (a + b + c) * b = 105
    zpie_element e_addsmul;
    zpie_init(&e_addsmul);
    zpie_element arr[3];
    arr[0] = a;
    arr[1] = b;
    arr[2] = c;
    int sz = 3;
    zpie_addsmul(&e_addsmul, &sz, arr, &b);

    // zpie_assert_equal: e_mul should equal a*b = 21
    zpie_element e_check;
    zpie_init(&e_check);
    zpie_mul(&e_check, &a, &b);
    zpie_assert_equal(&e_check, &e_mul);
}

Test(zpie, all_constraints)
{
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_all_constraints);
    zpie_proof p;
    zpie_generate_proof(&p, &test_all_constraints, &keys.pk);
    cr_assert(zpie_verify_proof(&test_all_constraints, &p, &keys.vk));
}

Test(zpie, serialization_roundtrip)
{
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_mimc_hash);
    zpie_store_setup(&keys);

    zpie_setup_keys keys2;
    zpie_read_setup(&keys2, &test_mimc_hash);

    // compare deserialized keys against originals
    cr_assert_eq(keys.pk.Ne, keys2.pk.Ne);
    cr_assert_eq(keys.pk.qap_size, keys2.pk.qap_size);
    cr_assert(mclBnGT_isEqual(&keys.vk.alphabetaT, &keys2.vk.alphabetaT));
    cr_assert(mclBnG2_isEqual(&keys.vk.gamma2, &keys2.vk.gamma2));
    cr_assert(mclBnG2_isEqual(&keys.vk.delta2, &keys2.vk.delta2));
    for (int i = 0; i < (nPublic + nConst); i++)
        cr_assert(mclBnG1_isEqual(&keys.vk.vk1[i], &keys2.vk.vk1[i]));
    cr_assert(mclBnG1_isEqual(&keys.pk.alpha1, &keys2.pk.alpha1));
    cr_assert(mclBnG2_isEqual(&keys.pk.beta2, &keys2.pk.beta2));
    cr_assert(mclBnG1_isEqual(&keys.pk.delta1, &keys2.pk.delta1));

    // prove with deserialized keys and verify
    zpie_proof p;
    zpie_generate_proof(&p, &test_mimc_hash, &keys2.pk);
    cr_assert(zpie_verify_proof(&test_mimc_hash, &p, &keys2.vk));
}

Test(zpie, invalid_proof_rejected)
{
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_single_constraint);
    zpie_proof p;
    zpie_generate_proof(&p, &test_single_constraint, &keys.pk);

    // zpie_proof should verify normally
    cr_assert(zpie_verify_proof(&test_single_constraint, &p, &keys.vk));

    // tamper with piA — add generator to it
    mclBnG1 g;
    mclBnG1_setStr(&g, GGEN, strlen(GGEN), 10);
    mclBnG1_add(&p.piA, &p.piA, &g);

    // tampered zpie_proof must not verify
    cr_assert_not(zpie_verify_proof(&test_single_constraint, &p, &keys.vk));
}

Test(zpie, cross_circuit_rejected)
{
    zpie_setup_keys keys_sc;
    zpie_perform_setup(&keys_sc, &test_single_constraint);
    zpie_setup_keys keys_mh;
    zpie_perform_setup(&keys_mh, &test_mimc_hash);

    zpie_proof p_sc;
    zpie_generate_proof(&p_sc, &test_single_constraint, &keys_sc.pk);

    // zpie_proof for single_constraint must not verify against mimc_hash's vk
    cr_assert_not(zpie_verify_proof(&test_single_constraint, &p_sc, &keys_mh.vk));
}

Test(zpie, deterministic_mimc)
{
    test_no_rand = 1;
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_mimc_hash);
    zpie_proof p;
    zpie_generate_proof(&p, &test_mimc_hash, &keys.pk);

    // verify the deterministic zpie_proof
    cr_assert(zpie_verify_proof(&test_mimc_hash, &p, &keys.vk));

    // generate a second time and check reproducibility
    zpie_proof p2;
    zpie_generate_proof(&p2, &test_mimc_hash, &keys.pk);
    cr_assert(mclBnG1_isEqual(&p.piA, &p2.piA));
    cr_assert(mclBnG2_isEqual(&p.piB2, &p2.piB2));
    cr_assert(mclBnG1_isEqual(&p.piC, &p2.piC));

    test_no_rand = 0;
}

Test(zpie, proof_serialization)
{
    zpie_setup_keys keys;
    zpie_perform_setup(&keys, &test_single_constraint);
    zpie_store_setup(&keys);

    zpie_proof p;
    zpie_generate_proof(&p, &test_single_constraint, &keys.pk);
    zpie_store_proof(&p);

    zpie_setup_keys keys2;
    zpie_read_setup(&keys2, &test_single_constraint);
    zpie_proof p2;
    zpie_read_proof(&p2);

    cr_assert(mclBnG1_isEqual(&p.piA, &p2.piA));
    cr_assert(mclBnG2_isEqual(&p.piB2, &p2.piB2));
    cr_assert(mclBnG1_isEqual(&p.piC, &p2.piC));

    cr_assert(zpie_verify_proof(&test_single_constraint, &p2, &keys2.vk));
}
