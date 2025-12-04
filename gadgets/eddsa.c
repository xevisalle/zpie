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
    element out[4];

    for (int i = 0; i < 4; ++i)
    {
        init(&out[i]);
    }

    element outPrivate[3];
    init_array(outPrivate, 3);

    element Bx, By;
    init(&Bx);
    init(&By);

    input(&Bx, B.x);
    input(&By, B.y);

    int arraySize = 5;
    element ram[arraySize];
    init_array(ram, arraySize);

    input(&ram[0], edsig.R.x);
    input(&ram[1], edsig.R.y);
    input(&ram[2], A.x);
    input(&ram[3], A.y);
    input(&ram[4], msg);

    element signature;
    init(&signature);
    input(&signature, edsig.S);

    multi_hash(outPrivate[0], ram, 5);

    int size = 254;
    element hBits[size];
    element sBits[size];

    init_array(hBits, size);
    init_array(sBits, size);

    to_bits(hBits, outPrivate[0], size);
    to_bits(sBits, signature, size - 1);

    mul_scalar(outPrivate[1], outPrivate[2], ram[2], ram[3], hBits, size);
    add_points(out[0], out[1], outPrivate[1], outPrivate[2], ram[0], ram[1]);

    mul_scalar(out[2], out[3], Bx, By, sBits, size - 1);

    assert_equal(&out[2], &out[0]);
    assert_equal(&out[3], &out[1]);
}

#endif