#include "../src/zpie.h"

int main()
{
    // we init the bulletproofs module, for 2 aggregated proofs of 64 bits
    bulletproof_init(64, 2);

    // we set some values to prove knowledge of, and compute the proof (../data/bulletproof.params)
    unsigned char *si[] = {"1234", "5678"};
    bulletproof_prove(si);

    // we verify the bulletproof (../data/bulletproof.params)
    if(bulletproof_verify()) printf("Bulletproof verified.\n");
    else printf("Bulletproof cannot be verified.\n");
}