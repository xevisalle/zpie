# ZPiE: Zero-knowledge Proofs in Embedded systems

A portable and efficient C library for developing Zero-Knowledge applications for embedded systems. 

**DISCLAIMER**: this library is currently **unstable**. Furthermore, **it has not gone through an exhaustive security analysis**, so it is not intended to be used in a production environment, only for academic purposes.

An academic paper about ZPiE has been published in the special issue *Recent Advances in Security, Privacy, and Applied Cryptography* of the journal *Mathematics* (2021) and can be found [here](https://doi.org/10.3390/math9202569).

## Overview

ZPiE supports the following Zero-Knowledge schemes, defined over the elliptic curves BN128 and BLS12-381:

- zk-SNARKs for arithmetic circuits. We support the [Groth'16](https://eprint.iacr.org/2016/260.pdf) scheme. ZPiE includes the following arithmetic circuits:
    - [EdDSA](https://eprint.iacr.org/2015/677.pdf) signature algorithm over [Baby JubJub](https://iden3-docs.readthedocs.io/en/latest/_downloads/33717d75ab84e11313cc0d8a090b636f/Baby-Jubjub.pdf) elliptic curve and BN128.
    - [MiMC-7](https://eprint.iacr.org/2016/492.pdf) hashing function (BN128 order).
- [Bulletproofs](https://eprint.iacr.org/2017/1066.pdf). We support range proofs (and aggregated range proofs).

In order to compute the circuit inputs for the above described circuits, you can use this [repository](https://github.com/xevisalle/cryptoolz).


## Install dependencies
ZPiE needs [GMP](https://gmplib.org/) and [MCL](https://github.com/herumi/mcl). To install them, and some other required dependencies, simply run:

```
sudo apt install libgmp-dev libcunit1-dev
git clone https://github.com/herumi/mcl
cd mcl
make -j4
```

## Test
ZPiE can be tested as follows:

```
git clone https://github.com/xevisalle/zpie
cd zpie
make test
```

## Benchmarks

You can compile a benchmarking application by running the following command, which will also prompt you with a set of available options to benchmark all the features of ZPiE:

```
make bench
```

## Compiling options

We can specify the elliptic curve to be used:

```
make bench CURVE=[OPTION]

Where [OPTION] can be:
BN128 (default)
BLS12_381
```

We can specify to run the code in multi-thread mode:

```
make bench MULTI=on
```

## zk-SNARKs for arithmetic circuits

Here there is an example on how to use zk-SNARKs. Copy the following snippet into a file (e.g. called `/src/main.c`):

```c
#include "zpie.h"

void circuit()
{
    // set a public element
    element out;
    init_public(&out);

    // set private elements
    element a, b;
    init(&a);
    init(&b);

    // input a value for the elements
    input(&a, "1234");
    input(&b, "5678");

    // apply a constraint multiplying such elements
    mul(&out, &a, &b);
}

int main()
{
    // we perform the setup
    setup_keys keys = perform_setup(&circuit); 

    // we generate a proof
    proof p = generate_proof(&circuit, keys.pk);

    // we verify the proof 
    if (verify_proof(&circuit, p, keys.vk)) 
        printf("Proof verified.\n");
    else 
        printf("Proof cannot be verified.\n");
}
```

And compile and execute using:

```
make MAIN=main && ./zpie
```

More circuit examples can be found in the `/src/tests.c` file.

## Bulletproofs

TBC.

## Cross-compile

### Build for x86_64

First we need to download and compile GMP for i386 64-bits:

```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar -xf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --enable-cxx ABI=64
make -j12
sudo make install
```

Then, we have to build MCL for i386 64-bits:

```
git clone https://github.com/herumi/mcl
cd mcl
make -j12 ARCH=x86_64
```

We finally build ZPiE:

```
make bench ARCH=x86_64
```

### Build for x86

First we need to download and compile GMP for i386 32-bits:

```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar -xf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --enable-cxx ABI=32
make -j12
sudo make install
```

Then, we have to build MCL for i386 32-bits:

```
git clone https://github.com/herumi/mcl
cd mcl
make -j12 ARCH=x86
```

We finally build ZPiE:

```
make bench ARCH=x86
```

### Build for ARM 64-bits

First we need to download and compile GMP for ARM 64-bits:

```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar -xf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --host=aarch64 ABI=64 CC=aarch64-linux-gnu-gcc
make -j12
sudo make install
```

Then, we have to build MCL for 32-bits:

```
git clone https://github.com/herumi/mcl
cd mcl
make -j12 CXX=aarch64-linux-gnu-g++ ARCH=aarch64 MCL_USE_GMP=0
```

We finally build ZPiE:

```
make bench ARCH=aarch64
```

### Build for ARM 32-bits

First we need to download and compile GMP for ARM 32-bits:

```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar -xf gmp-6.2.0.tar.xz
cd gmp-6.2.0
./configure --enable-cxx --host=armv6l ABI=32 CC=arm-linux-gnueabihf-gcc
make -j12
sudo make install
```

Then, we have to build MCL for ARM 32-bits:

```
git clone https://github.com/herumi/mcl
cd mcl
make -j12 CXX=arm-linux-gnueabihf-g++ ARCH=armv6l CFLAGS_USER="-I /usr/local/include" LDFLAGS="/usr/local/lib/libgmp.a /usr/local/lib/libgmpxx.a"
```

We finally build ZPiE:

```
make bench ARCH=arm
```
