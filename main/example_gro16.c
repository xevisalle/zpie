#include "../src/zpie.h"

int main()
{
	// we perform the setup (../data/provingkey.params and ../data/verifyingkey.params)
    init_setup();
	perform_setup();   

	// we generate a proof (../data/proof.params)
    init_prover();
    generate_proof();

    // we verify the proof (../data/proof.params)
    init_verifier();
    if (verify_proof()) printf("Proof verified.\n");
    else printf("Proof cannot be verified.\n");
}