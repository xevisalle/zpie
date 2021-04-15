#include "../src/zpie.h"

int main(int argc, char *argv[])
{   
    bench = 1; 
    if (argc < 2)
    {
        printf("******************* ZPiE v0.1 *******************\n");
        printf("USAGE: ./zpie [ACTIONS] [OPTIONS]\n\n");
        printf("[ACTIONS]:\n");
        printf("-s : Perform setup.\n");
        printf("-p : Generate proof.\n");
        printf("-v : Verify proof.\n\n");
        printf("[OPTIONS]\n");
        printf("-l : Activate operation logs.\n");

        exit(0);
    }

    if ((argc == 3) && (strcmp(argv[2], "-l") == 0)) logs = 1;

    init_setup();

    printf("******************* ZPiE v0.1 *******************\n");
    printf("--- Starting ZPiE...\n");
    printf("  |--- # of constraints: %d\n", N);
    printf("  |--- # of variables: %d\n", M);
    printf("  |--- # of public outputs: %d\n", nPublic);
    #ifdef MULTI
    printf("  |--- Multi-core execution: ON\n");
    #else
    printf("  |--- Multi-core execution: OFF\n");
    #endif

    #ifdef BN128
    printf("  |--- Elliptic curve: BN128\n");
    #else
    printf("  |--- Elliptic curve: BLS12_381\n");
    #endif

    if (strcmp(argv[1], "-s") == 0) perform_setup();
    else if (strcmp(argv[1], "-p") == 0)
    {
        init_prover();
        generate_proof();
    }
    else if (strcmp(argv[1], "-v") == 0)
    {
        init_verifier();
        verify_proof();
    }

    return 0;
}