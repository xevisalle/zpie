CC = gcc
OUT = build/libzpie.a
CAARCH64 = aarch64-linux-gnu-gcc
CARM = arm-linux-gnueabihf-gcc

.PHONY: zpie test bench clean
COMMON = -c src/zpie.c -o build/zpie.o -std=gnu17 -O3 -ffast-math -Wno-pointer-sign

MCLPATH = lib/mcl

MCLINCL = -I $(MCLPATH)/include
EXTLIB = -lcriterion -lm -lstdc++
MCLLIB = $(MCLPATH)/lib/lishe384_256.a $(MCLPATH)/lib/libmcl.a

EXTLIBMAC = -I /opt/homebrew/opt/libomp/include -I /opt/homebrew/include -L /opt/homebrew/lib -lcriterion -lstdc++

EXTLIBCROSS = -lstdc++

SRC = $(CURDIR)/src/*.c $(CURDIR)/gadgets/*.c

CURVE = BN128
ARCH = None
MULTI = off

ifeq ($(MULTI), on)
	MULTI_SET = -DMULTI_SET -fopenmp
	MCL_OMP = MCL_USE_OMP=1
endif 

zpie: $(SRC)
	mkdir -p build
	git submodule update --init
	cd lib/mcl && make -j16 $(MCL_OMP)
ifeq ($(ARCH), x86)
	$(CC) -m32 $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D$(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), x86_64)
	$(CC) -m64 $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D$(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), aarch64)
	$(CAARCH64) $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D$(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), arm)
	$(CARM) $(COMMON) $(MCLINCL) $(EXTLIBCROSS) -D$(CURVE) $(MULTI_SET)

else ifeq ($(shell uname), Darwin)
	$(CC) $(COMMON) $(MCLINCL) $(EXTLIBMAC) -D$(CURVE) $(MULTI_SET) -DIS_MAC_OS

else
	$(CC) $(COMMON) $(MCLINCL) $(EXTLIB) -D$(CURVE) $(MULTI_SET)

endif
	ar rc build/libzpie.a build/zpie.o
test: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	$(MAKE)
ifeq ($(shell uname), Darwin)
	$(CC) tests/test_zpie.c -o build/test_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIBMAC) -D$(CURVE) $(MULTI_SET) -DIS_MAC_OS
else
	$(CC) tests/test_zpie.c -o build/test_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D$(CURVE) $(MULTI_SET)
endif
	./build/test_zpie
bench: 
ifneq ("$(wildcard build/libzpie.a)","")
	rm build/libzpie.a
endif
	$(MAKE)
ifeq ($(shell uname), Darwin)
	$(CC) tests/bench_zpie.c -o build/bench_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIBMAC) -D$(CURVE) $(MULTI_SET) -DIS_MAC_OS
else
	$(CC) tests/bench_zpie.c -o build/bench_zpie build/libzpie.a $(MCLLIB) -I ./include $(MCLINCL) $(EXTLIB) -D$(CURVE) $(MULTI_SET)
endif
	./build/bench_zpie
clean:
	rm -rf data
	rm -rf build
	cd lib/mcl && make clean
