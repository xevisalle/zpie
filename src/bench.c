#include "zpie.h"

void circuit()
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

int main(int argc, char *argv[])
{   
    bench = 1; 
    if (argc < 2)
    {
        printf("******************* ZPiE v0.2 *******************\n");
        printf("USAGE: ./zpie [ACTIONS] [OPTIONS]\n\n");
        printf("[ACTIONS]:\n");
        printf("-s : Perform setup.\n");
        printf("-p : Generate proof.\n");
        printf("-v : Verify proof.\n");
        printf("-pbp <Nb> <Mc> : Generate bulletproof where Nb is the bit size and Mc the number of aggregated proofs.\n");
        printf("-vbp <Nb> <Mc> : Verify bulletproof where Nb is the bit size and Mc the number of aggregated proofs.\n\n");
        printf("[OPTIONS]\n");
        printf("-l : Activate operation logs.\n");

        exit(0);
    }

    if ((argc == 3) && (strcmp(argv[2], "-l") == 0)) logs = 1;

    printf("******************* ZPiE v0.2 *******************\n");

    if ((strcmp(argv[1], "-s") == 0) || (strcmp(argv[1], "-p") == 0) || (strcmp(argv[1], "-v") == 0))
    {
        init_setup();

        printf("--- Starting ZPiE - Groth'16...\n");
        printf("  |--- # of constraints: %d\n", N);
        printf("  |--- # of variables: %d\n", M);
        printf("  |--- # of public outputs: %d\n", nPublic);
    }
    else
    {
        printf("--- Starting ZPiE - Bulletproofs...\n");
        printf("  |--- # of bits : %s\n", argv[2]);
        printf("  |--- # of aggregated proofs: %s\n", argv[3]);       
    }
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


    if (strcmp(argv[1], "-s") == 0)
    {
        setupKeys keys = perform_setup();
        store_setup(keys);
    }
    else if (strcmp(argv[1], "-p") == 0)
    {
        setupKeys keys = read_setup();
        proof p = generate_proof(keys.pk);
        store_proof(p);
    }
    else if (strcmp(argv[1], "-v") == 0)
    {
        setupKeys keys = read_setup();
        proof p = read_proof();
        verify_proof(p, keys.vk);
    }
    else if (strcmp(argv[1], "-pbp") == 0)
    {
        int Nbu = atoi(argv[2]);
        int Mcu = atoi(argv[3]);

        bulletproof_init(Nbu, Mcu);
        unsigned char *si[Mcu];

        for (int i = 0; i < Mcu; i++)
        {
            si[i] = "2";
        }

        bulletproof_prove(si);
    }
    else if (strcmp(argv[1], "-vbp") == 0)
    {
        int Nbu = atoi(argv[2]);
        int Mcu = atoi(argv[3]);

        bulletproof_init(Nbu, Mcu);
        bulletproof_verify();
    }

    return 0;
}