FFTW=/opt/fftw/3.2.2/gnu
CC=gcc
CFLAGS=-O2 -Wall -g -I$(FFTW)/include 
LDFLAGS=-L$(FFTW)/lib -lfftw3_threads -lfftw3 -lm

all: becon

becon: becon.o io_grafic.o log.o
	$(CC) $(LDFLAGS) -o $@ $^

becon.o: becon.c becon.h
	$(CC) $(CFLAGS) -c $<

log.o: log.c log.h
	$(CC) $(CFLAGS) -c $<

io_grafic.o: io_grafic.c io_grafic.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o becon
