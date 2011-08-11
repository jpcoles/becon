#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "jpeglib.h"
#ifdef __cplusplus
}
#endif

#include "becon.h"
#include "frame_buffer.h"
#include "io_image.h"

int write_image(const char *fmt, int step, struct frame_buffer *fb)
{
    char *fname;

    asprintf(&fname, fmt, step);
    if (fname == NULL) return 1;
    sprintf(fname, fmt, step);

    FILE *fp = fopen(fname, "wb");
    if (fp != NULL)
    {
        write_frame_buffer(fp, fb);
        fclose(fp);
    }

    free(fname);

    return 0;
}

int capture_image(struct space *space, struct state *state, struct frame_buffer *fb)
{
    int32_t i,j,k;
    size_t idx = 0;

    double rho_max = 0;

    memset(fb->buf, 0, fb->width * fb->height * fb->samples);

    for (idx=0; idx < state->N; idx++)
        if (state->rho[idx][0] > rho_max) 
            rho_max = state->rho[idx][0];

    idx = 0;
    for (k=0; k < space->Nz; k++)
    for (j=0; j < space->Ny; j++)
    for (i=0; i < space->Nx; i++)
    {
        size_t ix = ((float)i)/space->Nx * fb->width;
        size_t iy = ((float)j)/space->Ny * fb->height;

        size_t o = fb->samples * (iy * fb->width + ix);

        int v = (state->rho[idx][0] / rho_max) * 255;
        //if (idx < space->Nx)
            //fprintf(stderr, "%ld %i %f %ld %ld\n", o, v, state->rho[idx][0], ix, iy);

        if (v > fb->buf[o+0])
        {
            fb->buf[o+0] = 
            fb->buf[o+1] = 
            fb->buf[o+2] = v;
            //fb->buf[o+3] = 1;
        }

        idx++;
    }

    return 0;
}
