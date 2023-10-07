#include "../src/zpie.h"

#include "../circuits/mimc.c"
#include "../circuits/eddsa.c"

// do x multiplications
void test1(int mulsize)
{
    element out;
    init_public(&out);

    element arr[mulsize];
    init_array(arr, mulsize);

    input(&arr[0], "213");
    input(&arr[1], "1122");

    for (int i = 2; i < mulsize; i++)
    {
        mul(&arr[i], &arr[1], &arr[i-1]);
    }

    mul(&out, &arr[0], &arr[mulsize-1]);
}

// verify an EdDSA signature
void test2()
{
    char *B1 = "5299619240641551281634865583518297030282874472190772894086521144482721001553";
    char *B2 = "16950150798460657717958625567821834550301663161624707787222815936182638968203";
    
    char *R1 = "1262948111445225057373438194818763405700457487429548371463214326190311895864"; 
    char *R2 = "12533500305127747239777484416561675628195562065959201739446841668623540883587";
    char *A1 = "21629779320182474195265732521833299809982444552305142529409236301104997786342";
    char *A2 = "9011812445381030664142622066218331845140881847034934166630871421746105699091";
    char *msg = "1234";

    char *signature = "2674591880888862378688383832785447197125897205360861957116147165712709455207";

    verify_eddsa(B1, B2, R1, R2, A1, A2, msg, signature);
}

// compute a MiMC hash
void test3()
{
    element h, x_in, k;

    init_public(&h);
    init(&x_in);
    init(&k);
    
    input(&x_in, "1234");
    input(&k, "112233445566");

    mimc7(&h, &x_in, &k);
}

void circuit()
{
    element out;
    init_public(&out);

    element arr0, arr1;
    init(&arr0);
    init(&arr1);

    input(&arr0, "213");
    input(&arr1, "1122");

    mul(&out, &arr0, &arr1);
}

int main()
{
	// we perform the setup
	setupKeys keys = perform_setup();  

    // we store and read the setup (../data/provingkey.params and ../data/verifyingkey.params)
    store_setup(keys);
    keys = read_setup();

	// we generate a proof
    proof p = generate_proof(keys.pk);

    // we verify the proof
    if (verify_proof(p, keys.vk)) printf("Proof verified.\n");
    else printf("Proof cannot be verified.\n");
}