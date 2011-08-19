#define _GNU_SOURCE
#include <stdlib.h>
#include <stdarg.h>
#include "log.h"
#include "io_tipsy.h"

int write_tipsy_grid(struct env *env, const char *fmt, ...)
{
    size_t idx;
    int32_t i,j,k;
    char *fname;
    va_list ap;

    struct tipsy_header hdr;
    struct tipsy_dm dm;

    va_start(ap, fmt);
    vasprintf(&fname, fmt, ap);
    va_end(ap);

    if (fname == NULL) return 1;

    va_start(ap, fmt);
    vsprintf(fname, fmt, ap);
    va_end(ap);

    FILE *fp = fopen(fname, "wb");
    if (fp != NULL)
    {
        Log("Writing %s...\n", fname);

        hdr.time = 0;
        hdr.nBodies = env->state.N;
        hdr.nDims = 3;
        hdr.nDark = hdr.nBodies;
        hdr.nStar = 0;
        hdr.nGas = 0;

        fwrite(&hdr, sizeof(hdr), 1, fp);

        for (k=0; k < env->space.Nx; k++)
        for (j=0; j < env->space.Ny; j++)
        for (i=0; i < env->space.Nz; i++)
        {
            dm.pos[0] = env->space.xmin + i*env->space.dx;
            dm.pos[1] = env->space.ymin + j*env->space.dx;
            dm.pos[2] = env->space.zmin + k*env->space.dx;

            dm.vel[0] = 
            dm.vel[1] = 
            dm.vel[2] = 0;

            dm.eps = 
            dm.phi = 
            dm.mass = 0;

            fwrite(&dm, sizeof(dm), 1, fp);
        }


        fclose(fp);
    }

    free(fname);

    return 0;
}
