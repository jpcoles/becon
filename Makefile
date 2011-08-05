FFTW=/opt/fftw/3.2.2/gnu
CC=gcc
CFLAGS=-O2 -Wall -g -I$(FFTW)/include 
LDFLAGS=-L$(FFTW)/lib -lfftw3_threads -lfftw3 -lm

all: becon

becon: becon.o
	$(CC) $(LDFLAGS) -o $@ $?

becon.o: becon.c becon.h
	$(CC) $(CFLAGS) -c $?

clean: 
	rm -f *.o becon
