FFTW=/opt/fftw/3.2.2/gnu
CC=gcc
CFLAGS=-O2 -Wall -g -I$(FFTW)/include 
CFLAGS+=-funroll-loops -ftree-vectorize -fno-omit-frame-pointer -fprefetch-loop-arrays -mssse3
LDFLAGS=-L$(FFTW)/lib -lfftw3_threads -lfftw3 -lm -ljpeg

all: becon

becon: becon.o io_grafic.o log.o io_image.o frame_buffer.o
	$(CC) $(LDFLAGS) -o $@ $^

becon.o: becon.c becon.h
	$(CC) $(CFLAGS) -c $<

log.o: log.c log.h
	$(CC) $(CFLAGS) -c $<

io_grafic.o: io_grafic.c io_grafic.h
	$(CC) $(CFLAGS) -c $<

io_image.o: io_image.c io_image.h
	$(CC) $(CFLAGS) -c $<

frame_buffer.o: frame_buffer.c frame_buffer.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o becon
