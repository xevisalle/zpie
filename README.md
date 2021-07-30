# ZPiE: Zero-knowledge Proofs in Embedded systems

A portable and efficient C library for developing Zero-Knowledge applications for embedded systems. 

**DISCLAIMER**: this library is currently **unstable**. Furthermore, **it has not gone through an exhaustive security analysis**, so it is not intended to be used in a production environment, only for academic purposes.

## Overview

ZPiE supports the following Zero-Knowledge schemes, defined over the elliptic curves BN128 and BLS12-381:

- zk-SNARKs for arithmetic circuits. We support the [Groth'16](https://eprint.iacr.org/2016/260.pdf) scheme. ZPiE includes the following arithmetic circuits:
    - [EdDSA](https://eprint.iacr.org/2015/677.pdf) signature algorithm over [Baby JubJub](https://iden3-docs.readthedocs.io/en/latest/_downloads/33717d75ab84e11313cc0d8a090b636f/Baby-Jubjub.pdf) elliptic curve and BN128.
    - [MiMC-7](https://eprint.iacr.org/2016/492.pdf) hashing function (BN128 order).
- [Bulletproofs](https://eprint.iacr.org/2017/1066.pdf). We support range proofs (and aggregated range proofs).

In order to compute the circuit inputs for the above described circuits, you can use this [repository](https://github.com/xevisalle/cryptoolz).


## Install dependencies
ZPiE needs [GMP](https://gmplib.org/) and [MCL](https://github.com/herumi/mcl). To install them, simply run:

```
sudo apt install libgmp-dev
git clone https://github.com/herumi/mcl
cd mcl
make -j4
```

## Compile
By default, a benchmarking application is compiled and executed as follows:

```
git clone https://github.com/xevisalle/zpie
cd zpie
make
./zpie
```

The last command will prompt you with a set of available options to benchmark all the features of ZPiE. Otherwise, you can compile your program as follows (which must be placed into `/main`):

```
make MAIN=example
```

### Compiling options

Optionally, following compiling targets can be specified (if not specified, 'single' target will be used):

```
make [target]

Where [target] can be:
single: single-threaded
multi: multi-threaded
```

We can also specify the multi-exponentiation algorithm to be used (by default 'bos-coster' is used):

```
make MULEXP=[OPTION]

Where [OPTION] can be:
BOSCOSTER_MULEXP: the Bos-Coster algorithm will be used.
NAIVE_MULEXP: serial multi-exponentiation will be performed.
MCL_MULEXP: the multi-exponentiation algorithm provided by the MCL library will be used.
```

We can also specify the elliptic curve to be used:

```
make CURVE=[OPTION]

Where [OPTION] can be:
BN128 (default)
BLS12_381
```


## zk-SNARKs for arithmetic circuits

### Circuit

Circuits in ZPiE are .c files, integrated into the binary. You can edit the file `circuit.c` to do modifications to the current circuit. A tutorial with all circuit options will be provided soon. Meanwhile, examples are provided into the `/circuits` folder.

### Setup, Prove, and Verify

Here there is an example on how to use zk-SNARKs (`/main/example_gro16`):

```c
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
```

## Bulletproofs

Here there is an example on how to use Bulletproofs (`/main/example_bp`):

```c
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
```

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
make ARCH=x86_64
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
make ARCH=x86
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
make ARCH=aarch64
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
make ARCH=arm
```
