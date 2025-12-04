CC = gcc
OUT = build/libzpie.a
CAARCH64 = aarch64-linux-gnu-gcc
CARM = arm-linux-gnueabihf-gcc
COMMON = -c src/zpie.c -o build/zpie.o -std=gnu99 -Ofast -Wno-unused-result -Wno-pointer-sign

MCLPATH = lib/mcl
GMPPATH = /usr/local

MCLINCL = -I $(MCLPATH)/include 
EXTLIB = -lgmp -lcriterion -lm -lstdc++
MCLLIB = $(MCLPATH)/lib/lishe384_256.a $(MCLPATH)/lib/libmcl.a

EXTLIBMAC = /opt/homebrew/lib/libgmp.a /opt/homebrew/opt/libomp/lib/libomp.a /opt/homebrew/opt/criterion/lib/libcriterion.a -I /opt/homebrew/opt/libomp/include -I /opt/homebrew/include -lm -lstdc++

EXTLIBCROSS = $(GMPPATH)/lib/libgmp.a -I $(GMPPATH)/include -lstdc++

SRC = $(shell pwd)/src/*.c $(shell pwd)/gadgets/*.c $(shell pwd)/include/*.h

CURVE = BN128
ARCH = None
MULTI = off

ifeq ($(MULTI), on)
	MULTI_SET = -D MULTI_SET -fopenmp
endif 

zpie: $(SRC)
	mkdir -p build
	git submodule update --init
	cd lib/mcl && make -j16
ifeq ($(ARCH), x86)
	$(CC) -m32 $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), x86_64)
	$(CC) -m64 $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), aarch64)
	$(CAARCH64) $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), arm)
	$(CARM) $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(shell uname), Darwin)
	$(CC) $(COMMON) $(MCLINCL) $(EXTLIBMAC) -D $(CURVE) $(MULTI_SET) IS_MAC_OS

else
	$(CC) $(COMMON) $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)

endif
	ar rcs build/libzpie.a build/zpie.o
test: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	make
	$(CC) tests/test_zpie.c -o build/test_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)
	./build/test_zpie
bench: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	make
	$(CC) tests/bench_zpie.c -o build/bench_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)
	./build/bench_zpie
clean:
	rm -rf data
	rm -rf build
