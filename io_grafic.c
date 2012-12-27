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

#define UNIT_LENGTH (56.36041080) // [m]
#define PC  (3.08568025e16) // [m]
#define KPC (1e3 * PC)
#define MPC (1e3 * KPC)

int alloc_grafic(const char *dirname, struct env *env)
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

    fclose(fp);

    env->cosmo.a_start = hdr.astart;
    env->cosmo.omega_m = hdr.omegam;
    env->cosmo.omega_v = hdr.omegav;
    env->cosmo.omega_r = 0;
    env->cosmo.omega_k = 1 - (env->cosmo.omega_m + env->cosmo.omega_v + env->cosmo.omega_r);
    env->cosmo.H0      = hdr.h0;
    env->cosmo.h       = hdr.h0 / 100.;

    double H2 = pow(env->cosmo.H0, 2); // * (env->cosmo.omega_v + ((((env->cosmo.omega_r / env->cosmo.a_start) + env->cosmo.omega_m) / env->cosmo.a_start + env->cosmo.omega_k) / env->cosmo.a_start / env->cosmo.a_start));
    env->cosmo.rho_crit = 3*H2 / (8*M_PI*env->consts.si.G) * (env->bec.m / (env->cosmo.H0 * env->consts.si.hbar));
    //env->cosmo.rho_crit /= 100000;

    env->consts.si.H0 = hdr.h0 / 3.08568025e19;

    //double Delta = hdr.dx / env->cosmo.h * 3.08568025e22;
    double Delta = hdr.dx * 3.08568025e22;
    env->eta = env->bec.m * pow(Delta,2) * pow(env->consts.si.H0,2) / env->consts.si.hbar;

#if 0
    fprintf(stderr, "dx %e\n", hdr.dx);
    fprintf(stderr, "h %e\n", env->cosmo.h);
    fprintf(stderr, "bec m %e\n", env->bec.m);
    fprintf(stderr, "H0 %e\n", env->consts.si.H0);
    fprintf(stderr, "Box %e\n", Delta * hdr.np1);
    fprintf(stderr, "Delta %e\n", Delta);
    fprintf(stderr, "ETA %f\n", env->eta);
    exit(0);
#endif


    env->space.Nx = hdr.np1;
    env->space.Ny = hdr.np2;
    env->space.Nz = hdr.np3;
    env->space.Nmax = MAX(env->space.Nx, MAX(env->space.Ny, env->space.Nz));

    env->space.dx = hdr.dx;
    env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/env->space.Nmax : 0;
    //env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nmax) : 0;
    //env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Ny) : 0;
    //env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nz) : 0;

    env->state.N = hdr.np1 * hdr.np2 * hdr.np3;

    env->state.psi = (fftw_complex*) fftw_malloc(sizeof(*env->state.psi) * env->state.N);
    env->state.phi = (fftw_complex*) fftw_malloc(sizeof(*env->state.phi) * env->state.N);
    env->state.rho = (fftw_complex*) fftw_malloc(sizeof(*env->state.rho) * env->state.N);

    Log("    Dimensions    %i x %i x %i\n", hdr.np1, hdr.np2, hdr.np3);
    Log("    Grid size     %f [comoving Mpc]\n", hdr.dx);
    Log("    Grid offsets  %f %f %f [comoving Mpc]\n", hdr.x1o, hdr.x2o, hdr.x3o);
    Log("    a start       %f\n", hdr.astart);
    Log("    omegam        %f\n", hdr.omegam);
    Log("    omegav        %f\n", hdr.omegav);
    Log("    omegar        %f\n", env->cosmo.omega_r);
    Log("    omegak        %f\n", env->cosmo.omega_k);
    Log("    rho_crit      %f\n", env->cosmo.rho_crit);
    Log("    H0            %f [km/s/Mpc]\n", hdr.h0);

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

int read_grafic(const char *dirname, struct env *env)
{
    size_t i;
    int retcode = 0;
    double rho_max=-DBL_MAX;
    double rho_min=DBL_MAX;

    float *arr = malloc(sizeof(*arr) * env->state.N);
    if (arr == NULL) 
    {
        retcode = 1;
        goto cleanup;
    }

    if (read_grafic_array(dirname, "/ic_deltab", arr, env->space.Nx, env->space.Ny, env->space.Nz))
    {
        retcode = 2;
        goto cleanup;
    }

    for (i=0; i < env->state.N; i++)
    {
        arr[i] += 1.0;
#if 0
        arr[i] *= env->cosmo.rho_crit;

        if (arr[i] < 0
#ifdef isnormal
		 || (arr[i] > 0 && !isnormal(arr[i]))
#endif
		)
        {
            Log("Bad input value %g\n", arr[i]);
            retcode = 3;
            goto cleanup;
        }
#endif
        env->state.phi[i][0] = 0;
        env->state.phi[i][1] = 0;

        env->state.psi[i][0] = sqrt(arr[i]);
        env->state.psi[i][1] = 0;

        env->state.rho[i][0] = arr[i];
        env->state.rho[i][1] = 0;

        if (arr[i] > rho_max) rho_max = arr[i];
        if (arr[i] < rho_min) rho_min = arr[i];
    }

    Log("rho_min is %f\n", rho_min);
    Log("rho_max is %f\n", rho_max);
    Log("1/sqrt(rho_max) is %f\n", 1/sqrt(rho_max));


cleanup:

    free(arr);

    return retcode;
}

