#include "../include/zpie.h"

int mulsize;

void bench_circuit()
{
    element a;
    init_public(&a);

    element arr[mulsize + 1];
    init_array(arr, mulsize + 1);

    input(&a, "1234");
    input(&arr[0], "5678");

    for (int i = 1; i <= mulsize; i++)
    {
        mul(&arr[i], &a, &arr[i - 1]);
    }
}

int main(int argc, char* argv[])
{
    bench = 1;
    if (argc < 3)
    {
        printf("******************* ZPiE v0.5 *******************\n");
        printf("USAGE: ./zpie [ACTIONS] [OPTIONS]\n\n");
        printf("[ACTIONS]:\n");
        printf("-s <c>: Perform setup of 'c' constraints.\n");
        printf("-p <c>: Generate proof of 'c' constraints.\n");
        printf("-v <c>: Verify proof of 'c' constraints.\n");
        printf("[OPTIONS]\n");
        printf("-l : Activate logs.\n");

        exit(0);
    }

    if ((argc == 4) && (strcmp(argv[3], "-l") == 0))
        logs = 1;

    printf("******************* ZPiE v0.5 *******************\n");

    if ((strcmp(argv[1], "-s") == 0) || (strcmp(argv[1], "-p") == 0) ||
        (strcmp(argv[1], "-v") == 0))
    {
        mulsize = strtol(argv[2], NULL, 10);
        init_setup(&bench_circuit);

        printf("--- Starting ZPiE - Groth'16...\n");
        printf("  |--- # of constraints: %d\n", N);
        printf("  |--- # of elements: %d\n", M);
        printf("  |--- # of public elements: %d\n", nPublic);
    }
    else
    {
        printf("--- Starting ZPiE - Bulletproofs...\n");
        printf("  |--- # of bits : %s\n", argv[2]);
        printf("  |--- # of aggregated proofs: %s\n", argv[3]);
    }
#ifdef MULTI_SET
    printf("  |--- Multi-core execution: ON\n");
#else
    printf("  |--- Multi-core execution: OFF\n");
#endif

#ifdef BN128
    printf("  |--- Elliptic curve: BN128\n");
#else
    printf("  |--- Elliptic curve: BLS12_381\n");
#endif

    if (strcmp(argv[1], "-s") == 0)
    {
        setup_keys keys = perform_setup(&bench_circuit);
        store_setup(&keys);
    }
    else if (strcmp(argv[1], "-p") == 0)
    {
        setup_keys keys = read_setup(&bench_circuit);
        proof p = generate_proof(&bench_circuit, &keys.pk);
        store_proof(&p);
    }
    else if (strcmp(argv[1], "-v") == 0)
    {
        setup_keys keys = read_setup(&bench_circuit);
        proof p = read_proof();
        verify_proof(&bench_circuit, &p, &keys.vk);
    }

    return 0;
}