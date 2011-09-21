#include <assert.h>
#include <math.h>
#include "becon.h"
#include "ic_spherical_collapse.h"
#include "log.h"

int ic_init_spherical_collapse(struct env *env)
{
    size_t idx, vidx;
    int32_t i,j,k;

    Log("Preparing spherical collapse ICs...\n");

    env->cosmo.a_start = 1;
    env->cosmo.H0 = 1;
    env->cosmo.omega_m = 1.0;
    env->cosmo.omega_v = 0.0;
    env->cosmo.omega_r = 0;
    env->cosmo.omega_k = 1 - (env->cosmo.omega_m + env->cosmo.omega_v + env->cosmo.omega_r);
    env->cosmo.rho_crit = 1;

    env->space.Nx = 128;
    env->space.Ny = 128;
    env->space.Nz = 1;
    env->space.Nmax = MAX(env->space.Nx, MAX(env->space.Ny, env->space.Nz));

    env->space.vNx = env->space.Nx;
    env->space.vNy = env->space.Ny;
    env->space.vNz = env->space.Nz;

    if (env->opts.with_vacuum)
    {
        if (env->space.vNx > 1) env->space.vNx *= 2;
        if (env->space.vNy > 1) env->space.vNy *= 2;
        if (env->space.vNz > 1) env->space.vNz *= 2;
    }
    
    env->space.dx = .1; //2*M_PI / env->space.Nmax / env->space.dk;
    env->space.dk = 2*M_PI / (env->space.dx * env->space.Nmax);

    assert(env->space.dx != 0);
    assert(env->space.dk != 0);

    env->space.xmin = -env->space.dx * env->space.Nx/2;
    env->space.xmax = -env->space.xmin;
    env->space.ymin = -env->space.dx * env->space.Ny/2;
    env->space.ymax = -env->space.ymin;
    env->space.zmin = -env->space.dx * env->space.Nz/2;
    env->space.zmax = -env->space.zmin;

    //env->space.dky = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Ny) : 0;
    //env->space.dkz = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nz) : 0;

    //fprintf(stderr, "%ld %i %i %i\n", idx, io, jo, ko);

    env->kmax2 = pow((env->space.Nx>>1)*env->space.dk,2)
               + pow((env->space.Ny>>1)*env->space.dk,2)
               + pow((env->space.Nz>>1)*env->space.dk,2);

    //Log("dt: %g\n", dt);
    Log("    dx: %g\n", env->space.dx);
    Log("    dk: %g\n", env->space.dk);
    Log(" kmax2: %g\n", env->kmax2);

    size_t N  = env->space.Nx * env->space.Ny * env->space.Nz;
    size_t vN = env->space.vNx * env->space.vNy * env->space.vNz;

    env->state.N  = N;
    env->state.vN = vN;

    Log("    %ld grid cells.\n", N);

    env->state.psi = (fftw_complex*) fftw_malloc(sizeof(*env->state.psi) * N);
    env->state.phi = (fftw_complex*) fftw_malloc(sizeof(*env->state.phi) * N);
    env->state.rho = (fftw_complex*) fftw_malloc(sizeof(*env->state.rho) * vN); /* vN *is* correct */

    if (env->opts.with_vacuum)
    {
        env->state.vphi = (fftw_complex*) fftw_malloc(sizeof(*env->state.vphi) * vN);
        env->state.S   = (fftw_complex*) fftw_malloc(sizeof(*env->state.S) * vN);
    }
    else
    {
        env->state.S   = NULL;
    }

    return 0;
}

int ic_spherical_collapse(struct env *env)
{
    int32_t i,j,k;
    size_t idx=0;

    struct space *sp = &env->space;

    double R = sp->Nx/4 * sp->dx;

    const double rho = 1e5;

    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->dx + sp->xmin;
        double y0 = j*sp->dx + sp->ymin;
        double z0 = k*sp->dx + sp->zmin;

        double r = sqrt(x0*x0 + y0*y0 + z0*z0);

        //if (0.8*R < r && r < R)
        if (r < R)
        {
            //env->state.rho[idx][0] = rho;
            //env->state.rho[idx][1] = 0;

            env->state.psi[idx][0] = sqrt(rho); // * pow(R-r,2));
            env->state.psi[idx][1] = 0;

            env->state.phi[idx][0] = 0;
            env->state.phi[idx][1] = 0;
        }
        else
        {
            //env->state.rho[idx][0] = 0;
            //env->state.rho[idx][1] = 0;

            env->state.psi[idx][0] = 0;
            env->state.psi[idx][1] = 0;

            env->state.phi[idx][0] = 0;
            env->state.phi[idx][1] = 0;
        }

        idx++;
    }

    env->consts.in.c = 0.1161270408e27;
    env->consts.in.H0 = .09032202799;
    env->consts.in.G  = 1;

    //env->consts.in.hbar = 0.00001290886701;
    //env->bec.m = 1e-4;

    env->consts.in.hbar = 10;
    env->bec.m = 1;

    env->drift_exp = env->consts.in.hbar / env->bec.m / 2;
    env->kick_exp  = env->bec.m / env->consts.in.hbar;

    //env->dt = 1e-5 * (1. / sqrt(env->consts.in.G * rho));
    env->dt = 1e-3 * 0.999 * M_PI_2 / env->kmax2 / env->drift_exp;

    fprintf(stderr, "%g should be < %g\n", 
        env->consts.in.hbar / env->bec.m,
        pow(sp->dx,2)*sqrt(env->consts.in.G * rho) / (4*M_PI));

    return 0;
}

