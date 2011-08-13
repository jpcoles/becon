#define _GNU_SOURCE
#include <assert.h>
#include <math.h>
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
#include "log.h"

int write_image(const char *fmt, int step, struct frame_buffer *fb)
{
    char *fname;

    asprintf(&fname, fmt, step);
    if (fname == NULL) return 1;
    sprintf(fname, fmt, step);

    FILE *fp = fopen(fname, "wb");
    if (fp != NULL)
    {
        Log("Writing %s...\n", fname);
        write_frame_buffer(fp, fb);
        fclose(fp);
    }

    free(fname);

    return 0;
}

int capture_image(struct space *space, struct state *state, struct frame_buffer *fb)
{
#if 0
    int32_t i,j,k;
    size_t idx = 0;

    double rho_max = 0;
    double rho_min = 1e20;

//#define Q(idx) (state->rho[idx][0])
#define Q(idx) (state->phi[idx][0])

    memset(fb->buf, 0, fb->width * fb->height * fb->samples);

    for (idx=0; idx < state->N; idx++)
    {
        double q = pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2);
        rho_max = fmax(q, rho_max);
        rho_min = fmin(q, rho_min);
    }

    idx = 0;
    #pragma omp parallel for private(j,i,idx)
    for (k=0; k < space->Nz; k++) { idx = k * space->Nx * space->Ny;
    for (j=0; j < space->Ny; j++)
    for (i=0; i < space->Nx; i++)
    {
        size_t ix = ((float)i)/space->Nx * fb->width;
        size_t iy = ((float)j)/space->Ny * fb->height;

        size_t o = fb->samples * (iy * fb->width + ix);

        double q = pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2);

        int v = ((q - rho_min) / fdim(rho_max, rho_min)) * 255;
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
    }

#endif
    return 0;
}

int capture_image_log(double min, double max, cmap_t cmap, struct space *space, struct state *state, struct frame_buffer *fb)
{
    int32_t i,j,k;
    size_t idx = 0;
    size_t idx2 = 0;

    memset(fb->buf, 0, fb->size * sizeof(*fb->buf));

    //#pragma omp parallel for private(j,i,idx)
    for (k=0; k < space->Nz; k++) { idx = k * space->Nx * space->Ny;
    for (j=0; j < space->Ny; j++)
    for (i=0; i < space->Nx; i++)
    {
        size_t ix = ((float)i)/space->Nx * fb->width;
        size_t iy = ((float)j)/space->Ny * fb->height;

        size_t o = iy * fb->width + ix;

        double q = pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2);
        double v = (log10(q) - log10(min)) / (log10(max) - log10(min));

        assert(v == 0 || isnormal(v));
        if (v < 0) v = 0;

        if (v > fb->buf[o])
            fb->buf[o] = v;

        idx++;
    }
    }

    #pragma omp parallel for private(idx2)
    for (idx=0; idx < fb->size; idx++)
    {
        int r,g,b;
        idx2 = idx * fb->samples;

        cmap(fb->buf[idx], &r, &g, &b);
        fb->img[idx2+0] = r;
        fb->img[idx2+1] = g;
        fb->img[idx2+2] = b;
    }

    return 0;
}
