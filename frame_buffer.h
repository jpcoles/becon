#ifndef FRAME_BUFFER_H
#define FRAME_BUFFER_H

#include <sys/types.h>
#include <stdio.h>
#include "jpeglib.h"

struct frame_buffer
{
    int32_t width, height;
    char samples;

    JSAMPROW buf;
};

int alloc_frame_buffer(struct frame_buffer *fb, int32_t width, int32_t height);
int free_frame_buffer(struct frame_buffer *fb);
void write_frame_buffer(FILE *outfile, struct frame_buffer *fb);

#endif
