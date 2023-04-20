CC=gcc
CFLAGS=-c -Ofast -funroll-loops -ffast-math -mavx2 -march=native -fPIC
LDFLAGS=-lm

all:
	$(CC) $(CFLAGS) kepler_scatteringPy.c -o solarsystemPy.o
	$(CC) -L. -shared solarsystemPy.o -o solarsystemPy.so -lm
	#$(CC) $(CFLAGS) kepler_hermite2body.c -o twobodyPy.o
	#$(CC) -L. -shared twobodyPy.o -o twobodyPy.so -lm


clean:
	rm -rf *.so *.o

