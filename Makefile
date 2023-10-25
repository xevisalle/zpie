MAIN = test
OUT = zpie
CC = gcc
CAARCH64 = aarch64-linux-gnu-gcc
CARM = arm-linux-gnueabihf-gcc
COMMON = src/$(MAIN).c -o $(OUT) -std=gnu99 -Ofast -Wno-unused-result
MCLPATH = ../mcl
GMPPATH = /usr/local
LIB = $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a -I $(MCLPATH)/include -lgmp -lcunit -lm -lstdc++
LIBMAC = /opt/homebrew/lib/libgmp.a /opt/homebrew/opt/libomp/lib/libomp.a /opt/homebrew/opt/cunit/lib/libcunit.a $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a -I /opt/homebrew/opt/libomp/include -I /opt/homebrew/include -I $(MCLPATH)/include -lm -lstdc++
LIBCROSS = $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a $(GMPPATH)/lib/libgmp.a -I $(MCLPATH)/include -I $(GMPPATH)/include -lstdc++
SRC = $(shell pwd)/src/*.c $(shell pwd)/circuits/*.c $(shell pwd)/src/*.h

MULEXP = MCL_MULEXP
CURVE = BN128
ARCH = None
MULTI = off

ifeq ($(MULTI), on)
	MULTI_SET = -D MULTI_SET -fopenmp
endif 

zpie: $(SRC)
ifeq ($(ARCH), x86)
	$(CC) -m32 $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), x86_64)
	$(CC) -m64 $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), aarch64)
	$(CAARCH64) $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

else ifeq ($(ARCH), arm)
	$(CARM) $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

else ifeq ($(shell uname), Darwin)
	$(CC) $(COMMON) $(LIBMAC) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

else
	$(CC) $(COMMON) $(LIB) -D $(MULEXP) -D $(CURVE) $(MULTI_SET)

endif
test: 
ifneq ("$(wildcard zpie)","")
	rm zpie
endif
	make MAIN=tests
	./zpie
bench: 
ifneq ("$(wildcard zpie)","")
	rm zpie
endif
	make MAIN=bench
	./zpie
clean:
	rm -rf data
ifneq ("$(wildcard zpie)","")
	rm zpie
endif