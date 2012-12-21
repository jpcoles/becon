FFTW_LDFLAGS=-L/usr/mpi/gcc/openmpi-1.4.3/lib64 -lfftw3 -lm 
FFTW_CFLAGS=-I/usr/mpi/gcc/openmpi-1.4.3/include

#FFTW_LDFLAGS+=-L/usr/local/lib #-lfftw3 -lm 
#FFTW_CFLAGS+=-I/usr/local/include
FFTW_LDFLAGS+=-L $(FFTW_LIB)
FFTW_CFLAGS+=-L $(FFTW_INCLUDE)

#FFTW_LDFLAGS=$(shell PKG_CONFIG_PATH=/opt/fftw/3.2.2/gnu/lib/pkgconfig pkg-config fftw3 --libs)
#FFTW_CFLAGS=$(shell PKG_CONFIG_PATH=/opt/fftw/3.2.2/gnu/lib/pkgconfig pkg-config fftw3 --cflags)

CFLAGS+=-L/usr/include/libpng
LDFLAGS+=-I/usr/lib -lpng

CC=gcc
CFLAGS+=-O3 -Wall -g -fopenmp $(FFTW_CFLAGS)
CFLAGS+=-funroll-loops -ftree-vectorize -fno-omit-frame-pointer -fprefetch-loop-arrays -mssse3
LDFLAGS+=-lfftw3_omp $(FFTW_LDFLAGS) -lgomp
#LDFLAGS=-lpthread -lfftw3_threads $(FFTW_LDFLAGS) -ljpeg -lgomp

#CFLAGS+=-I/usr/X11/include
#LDFLAGS+=-L/usr/X11/lib -lpng
#LDFLAGS+=-ljpeg

CFLAGS+=-I/opt/local/include
LDFLAGS+=-L/opt/local/lib -ljpeg

OBJS=becon.o io_grafic.o log.o io_image.o frame_buffer.o cmap.o io_arrays.o io_tipsy.o ic_spherical_collapse.o ic_infinite_sheet.o \
	analysis.o matrix.o

all: becon

becon: $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS) 

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o becon
