#include <math.h>
#include "ic_spherical_collapse.h"
#include "log.h"

int ic_init_spherical_collapse(struct env *env)
{
    Log("Preparing spherical collapse ICs...\n");

    env->cosmo.a_start = 1;
    env->cosmo.H0 = 1;
    env->cosmo.omega_m = 1.0;
    env->cosmo.omega_v = 0.0;
    env->cosmo.omega_r = 0;
    env->cosmo.omega_k = 1 - (env->cosmo.omega_m + env->cosmo.omega_v + env->cosmo.omega_r);
    env->cosmo.rho_crit = 1;

    env->space.Nx =
    env->space.Ny =
    env->space.Nz = 128;

    env->space.dx = 1. / env->space.Nz;

    env->space.xmin = -env->space.dx * env->space.Nx/2;
    env->space.xmax = -env->space.xmin;
    env->space.ymin = -env->space.dx * env->space.Ny/2;
    env->space.ymax = -env->space.ymin;
    env->space.zmin = -env->space.dx * env->space.Nz/2;
    env->space.zmax = -env->space.zmin;

    env->space.dkx = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nx) : 0;
    env->space.dky = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Ny) : 0;
    env->space.dkz = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nz) : 0;

    //fprintf(stderr, "%ld %i %i %i\n", idx, io, jo, ko);

    env->kmax2 = pow((env->space.Nx>>1)*env->space.dkx,2)
               + pow((env->space.Ny>>1)*env->space.dky,2)
               + pow((env->space.Nz>>1)*env->space.dkz,2);

    //Log("dt: %g\n", dt);
    Log("    dx: %g\n", env->space.dx);
    Log("    dk: %g %g %g\n", env->space.dkx, env->space.dky, env->space.dkz);

    size_t N = env->space.Nx * env->space.Ny * env->space.Nz;

    env->state.N = N;

    Log("    %ld grid cells.\n", N);

    env->state.psi = (fftw_complex*) fftw_malloc(sizeof(*env->state.psi) * N);
    env->state.phi = (fftw_complex*) fftw_malloc(sizeof(*env->state.phi) * N);
    env->state.rho = (fftw_complex*) fftw_malloc(sizeof(*env->state.rho) * N);
    return 0;
}

int ic_spherical_collapse(struct env *env)
{
    int32_t i,j,k;
    size_t idx=0;

    struct space *sp = &env->space;

    double R = sp->Nx/16 * sp->dx;

    const double rho = 1e9;

    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->dx + sp->xmin;
        double y0 = j*sp->dx + sp->ymin;
        double z0 = k*sp->dx + sp->zmin;

        double r = sqrt(x0*x0 + y0*y0 + z0*z0);

        if (r < R)
        {
            env->state.rho[idx][0] = rho;
            env->state.rho[idx][1] = 0;

            env->state.psi[idx][0] = sqrt(rho);
            env->state.psi[idx][1] = 0;

            env->state.phi[idx][0] = 0;
            env->state.phi[idx][1] = 0;
        }
        else
        {
            env->state.rho[idx][0] = 0;
            env->state.rho[idx][1] = 0;

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

    env->consts.in.hbar = 0.00001290886701;
    env->bec.m = 1e-4;

    env->drift_exp = env->consts.in.hbar / env->bec.m / 2;
    env->kick_exp  = env->bec.m / env->consts.in.hbar;

    env->dt = 1e-2 * (1. / sqrt(env->consts.in.G * rho));

    fprintf(stderr, "%g should be < %g\n", 
        env->consts.in.hbar / env->bec.m,
        pow(sp->dx,2)*sqrt(env->consts.in.G * rho) / (4*M_PI));

    return 0;
}

