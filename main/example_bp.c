// Bulletproofs are still experimental code

#include "../src/zpie.h"

int main()
{
    // we init the bulletproofs module
    bulletproof_init();

    // we set some values to prove knowledge of, and compute the proof (../data/bulletproof.params)
    unsigned char *si[] = {"1234", "5678"};
    bulletproof_prove(si);

    // we verify the bulletproof (../data/bulletproof.params)
    bulletproof_verify();   
}