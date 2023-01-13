CC=gcc
CFLAGS=-c -Ofast -mavx2 -march=native -fPIC
LDFLAGS=-lm

all:
	$(CC) $(CFLAGS) kepler_scatteringPy.c -o keplerPy.o
	$(CC) -L. -shared keplerPy.o -o keplerPy.so -lm
	$(CC) $(CFLAGS) kepler_hermite2body.c -o twobodyPy.o
	$(CC) -L. -shared twobodyPy.o -o twobodyPy.so -lm


clean:
	rm -rf keplerPy.so keplerPy.o

