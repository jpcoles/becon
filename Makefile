FFTW_LDFLAGS=$(shell PKG_CONFIG_PATH=/opt/fftw/3.2.2/gnu/lib/pkgconfig pkg-config fftw3 --libs)
FFTW_CFLAGS=$(shell PKG_CONFIG_PATH=/opt/fftw/3.2.2/gnu/lib/pkgconfig pkg-config fftw3 --cflags)
CC=gcc
CFLAGS=-O3 -Wall -g -fopenmp $(FFTW_CFLAGS)
CFLAGS+=-funroll-loops -ftree-vectorize -fno-omit-frame-pointer -fprefetch-loop-arrays -mssse3
LDFLAGS=-lpthread -lfftw3_threads $(FFTW_LDFLAGS) -ljpeg -lgomp

all: becon

becon: becon.o io_grafic.o log.o io_image.o frame_buffer.o cmap.o
	$(CC) -o $@ $^ $(LDFLAGS) 

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

cmap.o: cmap.c cmap.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o becon
