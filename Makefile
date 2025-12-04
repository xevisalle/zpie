CC = gcc
OUT = build/libzpie.a
CAARCH64 = aarch64-linux-gnu-gcc
CARM = arm-linux-gnueabihf-gcc
COMMON = -c src/zpie.c -o build/zpie.o -std=gnu99 -Ofast -Wno-unused-result -Wno-pointer-sign

MCLPATH = ../mcl
GMPPATH = /usr/local

MCLINCL = -I $(MCLPATH)/include 
EXTLIB = -lgmp -lcunit -lm -lstdc++
MCLLIB = $(MCLPATH)/lib/lishe384_256.a $(MCLPATH)/lib/libmcl.a

LIBMAC = /opt/homebrew/lib/libgmp.a /opt/homebrew/opt/libomp/lib/libomp.a /opt/homebrew/opt/cunit/lib/libcunit.a $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a -I /opt/homebrew/opt/libomp/include -I /opt/homebrew/include -I $(MCLPATH)/include -lm -lstdc++
LIBCROSS = $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a $(GMPPATH)/lib/libgmp.a -I $(MCLPATH)/include -I $(GMPPATH)/include -lstdc++
SRC = $(shell pwd)/src/*.c $(shell pwd)/circuits/*.c $(shell pwd)/include/*.h

CURVE = BN128
ARCH = None
MULTI = off

ifeq ($(MULTI), on)
	MULTI_SET = -D MULTI_SET -fopenmp
endif 

zpie: $(SRC)
	mkdir -p build
ifeq ($(ARCH), x86)
	$(CC) -m32 $(COMMON) $(LIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), x86_64)
	$(CC) -m64 $(COMMON) $(LIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), aarch64)
	$(CAARCH64) $(COMMON) $(LIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), arm)
	$(CARM) $(COMMON) $(LIBCROSS) -D $(CURVE) $(MULTI_SET)

else ifeq ($(shell uname), Darwin)
	$(CC) $(COMMON) $(LIBMAC) -D $(CURVE) $(MULTI_SET) -D IS_MAC_OS

else
	$(CC) $(COMMON) $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)
	ar rcs build/libzpie.a build/zpie.o

endif
test: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	make
	$(CC) src/tests.c -o build/zpie-tests build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)
	./build/zpie-tests
bench: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	make
	$(CC) src/bench.c -o build/zpie-bench build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D $(CURVE) $(MULTI_SET)
	./build/zpie-bench
clean:
	rm -rf data
ifneq ("$(wildcard zpie)","")
	rm zpie
endif