#include <float.h>
#include <math.h>
#include <sys/types.h>
#include <alloca.h>
#include <string.h>
#include <stdlib.h>

#include "becon.h"
#include "log.h"
#include "io_grafic.h"

struct header
{
    int32_t rec_begin;
    int32_t np1, np2, np3;
    float dx;
    float x1o, x2o, x3o;
    float astart, omegam, omegav;
    float h0;
    int32_t rec_end;
} hdr;

static int read_grafic_array(const char *dirname, const char *fname, float *a, size_t np1, size_t np2, size_t np3);

int alloc_grafic(const char *dirname, struct state *state, struct space *space)
{
    FILE *fp;
    char *fname = "/ic_deltab";
    char *fpath = alloca(strlen(dirname) + strlen(fname) + 1);

    sprintf(fpath, "%s%s", dirname, fname);

    Log("Reading %s\n", fpath);

    fp = fopen(fpath, "rb");
    if (fp == NULL) return 1;

    fread(&hdr, sizeof(hdr), 1, fp);

    if (hdr.rec_begin != hdr.rec_end)
    {
        Log("FORTRAN record lengths don't match in header.\n");
        return 1;
    }

    Log("    Dimensions    %i x %i x %i\n", hdr.np1, hdr.np2, hdr.np3);
    Log("    Grid size     %f [comoving Mpc]\n", hdr.dx);
    Log("    Grid offsets  %f %f %f [comoving Mpc]\n", hdr.x1o, hdr.x2o, hdr.x3o);
    Log("    a start       %f\n", hdr.astart);
    Log("    omegam        %f\n", hdr.omegam);
    Log("    omegav        %f\n", hdr.omegav);
    Log("    H0            %f [km/s/Mpc]\n", hdr.h0);

    fclose(fp);

    space->Nx = hdr.np1;
    space->Ny = hdr.np2;
    space->Nz = hdr.np3;

    space->drx =
    space->dry =
    space->drz = hdr.dx;

    space->dkx = fabs(space->drx) != 0 ? 2*M_PI/(space->drx * space->Nx) : 0;
    space->dky = fabs(space->dry) != 0 ? 2*M_PI/(space->dry * space->Ny) : 0;
    space->dkz = fabs(space->drz) != 0 ? 2*M_PI/(space->drz * space->Nz) : 0;

    state->N = hdr.np1 * hdr.np2 * hdr.np3;

    state->psi = (fftw_complex*) fftw_malloc(sizeof(*state->psi) * state->N);
    state->phi = (fftw_complex*) fftw_malloc(sizeof(*state->phi) * state->N);
    state->rho = (fftw_complex*) fftw_malloc(sizeof(*state->rho) * state->N);

    return 0;
}

static int read_grafic_array(const char *dirname, const char *fname, float *a, size_t np1, size_t np2, size_t np3)
{
    FILE *fp;
    //char *fname = "/ic_deltab";
    char *fpath = alloca(strlen(dirname) + strlen(fname) + 1);

    int retcode = 0;

    sprintf(fpath, "%s%s", dirname, fname);

    fp = fopen(fpath, "rb");
    if (fp == NULL) 
    {
        Log("Unable to open %s\n", fpath);
        retcode = 1;
        goto cleanup;
    }

    Log("Loading %s\n", fpath);

    fread(&hdr, sizeof(hdr), 1, fp);

    if (hdr.rec_begin != hdr.rec_end)
    {
        Log("    FORTRAN record lengths don't match in header.\n");
        retcode = 2;
        goto cleanup;
    }

    size_t i,j;
    int32_t rec_begin, rec_end;
    for (i=0; i < np3; i++)
    {
        fread(&rec_begin, sizeof(rec_begin), 1, fp);

        for (j=0; j < np2; j++)
        {
            if (fread(a, sizeof(*a), np1, fp) != np1)
            {
                Log("    Short read.");
                retcode = 3;
                goto cleanup;
            }

            if (ferror(fp))
            {
                Log("    File error.\n");
                retcode = 4;
                goto cleanup;
            }

            a += np1;
        }

        fread(&rec_end, sizeof(rec_end), 1, fp);

        if (rec_begin != rec_end)
        {
            Log("    FORTRAN record lengths don't match in header.\n");
            retcode = 2;
            goto cleanup;
        }
    }

cleanup:

    switch (retcode)
    {
        case 1:
            break;
        default:
            fclose(fp);
            break;
    }

    return retcode;
}

int read_grafic(const char *dirname, struct state *state, struct space *space)
{
    size_t i;
    int retcode = 0;
    double rho_max=-DBL_MAX;
    double rho_min=DBL_MAX;

    float *a = malloc(sizeof(*a) * state->N);
    if (a == NULL) 
    {
        retcode = 1;
        goto cleanup;
    }

    if (read_grafic_array(dirname, "/ic_deltab", a, space->Nx, space->Ny, space->Nz))
    {
        retcode = 2;
        goto cleanup;
    }

    for (i=0; i < state->N; i++)
    {
        a[i] += 1.0;

    //a[i] *= 100;

        if (a[i] < 0
#ifdef isnormal
		 || (a[i] > 0 && !isnormal(a[i]))
#endif
		)
        {
            Log("Bad input value %g\n", a[i]);
            retcode = 3;
            goto cleanup;
        }
        state->phi[i][0] = 0;
        state->phi[i][1] = 0;

        state->psi[i][0] = sqrt(a[i]);
        state->psi[i][1] = 0;

        state->rho[i][0] = a[i];
        state->rho[i][1] = 0;

        if (a[i] > rho_max) rho_max = a[i];
        if (a[i] < rho_min) rho_min = a[i];
    }

    Log("rho_min is %f\n", rho_min);
    Log("rho_max is %f\n", rho_max);
    Log("1/sqrt(rho_max) is %f\n", 1/sqrt(rho_max));


cleanup:

    free(a);

    return retcode;
}

