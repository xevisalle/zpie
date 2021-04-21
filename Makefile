MAIN = bench
OUT = zpie
CC = gcc
CAARCH64 = aarch64-linux-gnu-gcc
CARM = arm-linux-gnueabihf-gcc
COMMON = main/$(MAIN).c -o $(OUT) -std=gnu99 -Ofast
MCLPATH = ../mcl
GMPPATH = /usr/local
LIB = $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a -I $(MCLPATH)/include -lgmp -lstdc++
LIBCROSS = $(MCLPATH)/lib/libmclbn384_256.a $(MCLPATH)/lib/libmcl.a $(GMPPATH)/lib/libgmp.a -I $(MCLPATH)/include -I $(GMPPATH)/include -lstdc++
SRC = $(shell pwd)/src/*.c $(shell pwd)/circuits/*.c $(shell pwd)/main/*.c $(shell pwd)/src/*.h

MULEXP = BOSCOSTER_MULEXP
CURVE = BN128
ARCH = None

zpie: $(SRC)
ifeq ($(ARCH), x86)
	$(CC) -m32 $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI)

else ifeq ($(ARCH), x86_64)
	$(CC) -m64 $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI)

else ifeq ($(ARCH), aarch64)
	$(CAARCH64) $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI)

else ifeq ($(ARCH), arm)
	$(CARM) $(COMMON) $(LIBCROSS) -D $(MULEXP) -D $(CURVE) $(MULTI)

else
	$(CC) $(COMMON) $(LIB) -D $(MULEXP) -D $(CURVE) $(MULTI)

endif
	
single: 
ifneq ("$(wildcard zpie)","")
	rm zpie
endif
	make MULTI=""
multi:
ifneq ("$(wildcard zpie)","")
	rm zpie
endif
	make MULTI="-D MULTI -fopenmp"
clean:
	rm -rf data
ifneq ("$(wildcard zpie)","")
	rm zpie
endif