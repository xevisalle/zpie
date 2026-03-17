#ifndef EDDSA_C
#define EDDSA_C
#include "mimc.c"
#include "twisted_edwards.c"

typedef struct
{
    point R;
    char* S;
} eddsa_signature;

void verify_eddsa(eddsa_signature edsig, point B, point A, char* msg)
{
    zpie_element out[4];

    for (int i = 0; i < 4; ++i)
    {
        zpie_init(&out[i]);
    }

    zpie_element outPrivate[3];
    zpie_init_array(outPrivate, 3);

    zpie_element Bx, By;
    zpie_init(&Bx);
    zpie_init(&By);

    zpie_input(&Bx, B.x);
    zpie_input(&By, B.y);

    int arraySize = 5;
    zpie_element ram[arraySize];
    zpie_init_array(ram, arraySize);

    zpie_input(&ram[0], edsig.R.x);
    zpie_input(&ram[1], edsig.R.y);
    zpie_input(&ram[2], A.x);
    zpie_input(&ram[3], A.y);
    zpie_input(&ram[4], msg);

    zpie_element signature;
    zpie_init(&signature);
    zpie_input(&signature, edsig.S);

    multi_hash(outPrivate[0], ram, 5);

    int size = 254;
    zpie_element hBits[size];
    zpie_element sBits[size];

    zpie_init_array(hBits, size);
    zpie_init_array(sBits, size);

    to_bits(hBits, outPrivate[0], size);
    to_bits(sBits, signature, size - 1);

    mul_scalar(outPrivate[1], outPrivate[2], ram[2], ram[3], hBits, size);
    add_points(out[0], out[1], outPrivate[1], outPrivate[2], ram[0], ram[1]);

    mul_scalar(out[2], out[3], Bx, By, sBits, size - 1);

    zpie_assert_equal(&out[2], &out[0]);
    zpie_assert_equal(&out[3], &out[1]);
}

#endif