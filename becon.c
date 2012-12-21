#ifdef _OPENMP
#include <omp.h>
#endif

#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fftw3.h>

#include "becon.h"
#include "log.h"
#include "io_grafic.h"
#include "io_arrays.h"
#include "io_tipsy.h"
#include "io_image.h"
#include "frame_buffer.h"
#include "cmap.h"
#include "matrix.h"
#include "analysis.h"

#include "ic_spherical_collapse.h"
#include "ic_infinite_sheet.h"



//==============================================================================
//                                ic_init_test1D
//==============================================================================
void ic_init_test1D(struct env *env)
{
    Log("Preparing test1D initial conditions...\n");

    env->cosmo.a_start = 1;
    env->cosmo.H0 = 1;
    env->cosmo.omega_m = 0.279;
    env->cosmo.omega_v = 0.721;
    env->cosmo.omega_r = 0;
    env->cosmo.omega_k = 1 - (env->cosmo.omega_m + env->cosmo.omega_v + env->cosmo.omega_r);
    env->cosmo.rho_crit = 10;

    env->space.Nx = 1 << 10;
    env->space.Ny = 1;
    env->space.Nz = 1;

    env->space.Nmax = MAX(env->space.Nx, MAX(env->space.Ny, env->space.Nz));

    env->space.dx = 1; //(env->space.xmax-env->space.xmin) / env->space.Nx;

    env->space.xmin = -env->space.dx * env->space.Nx/2;
    env->space.xmax = -env->space.xmin;
    env->space.ymin = 0.0;
    env->space.ymax = 0.0;
    env->space.zmin = 0.0;
    env->space.zmax = 0.0;

    env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nmax) : 0;
    //env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Ny) : 0;
    //env->space.dk = fabs(env->space.dx) != 0 ? 2*M_PI/(env->space.dx * env->space.Nz) : 0;

    env->kmax2 = pow((env->space.Nx>>1)*env->space.dk,2)
               + pow((env->space.Ny>>1)*env->space.dk,2)
               + pow((env->space.Nz>>1)*env->space.dk,2);

    //Log("dt: %g\n", dt);
    Log("    kmax2: %g\n", env->kmax2);
    Log("    dx: %g\n", env->space.dx);
    Log("    dk: %g\n", env->space.dk);

    size_t N = env->space.Nx * env->space.Ny * env->space.Nz;

    env->state.N = N;

    Log("    %ld grid cells.\n", N);

    env->state.psi = (fftw_complex*) fftw_malloc(sizeof(*env->state.psi) * N);
    env->state.phi = (fftw_complex*) fftw_malloc(sizeof(*env->state.phi) * N);
    env->state.rho = (fftw_complex*) fftw_malloc(sizeof(*env->state.rho) * N);

}

//==============================================================================
//                                  ic_test1D
//==============================================================================
void ic_test1D(struct env *env)
{
    int32_t x;

    Log("Creating test1D initial conditions...\n");

    const double frac = 1./2.7;
    for (x=0; x < env->space.Nx; x++)
    {
        double rho = env->cosmo.rho_crit * (frac*env->space.Nx < x && x < (1-frac)*env->space.Nx);
        //double rho = env->cosmo.rho_crit * exp(-30 * pow(env->space.xmin + x*env->space.dx,2) / pow(env->space.xmax-env->space.xmin,2));

//      if (x < env->space.Nx/8 || x > env->space.Nx - (env->space.Nx / 8))
//          rho = 0;

        env->state.psi[x][0] = sqrt(rho);
        //fprintf(stderr, "%f\n", env->state.psi[x][0]);

        //s->psi[x][0] = cos(2*M_PI * x/sp->Nx);
        //s->psi[x][0] = 1e-5 * cos(2*M_PI * x/L);
        env->state.psi[x][1] = 0; //sin(2*M_PI * x/L);

        env->state.phi[x][0] = 0;
        env->state.phi[x][1] = 0;

        env->state.rho[x][0] = rho;
        env->state.rho[x][1] = 0;
    }
    env->consts.in.c = 0.1161270408e27;
    env->consts.in.H0 = .09032202799;
    env->consts.in.G  = 1;

    env->consts.in.hbar = 10;
    env->bec.m = 1;

    env->drift_exp = env->consts.in.hbar / env->bec.m / 2;
    env->kick_exp  = env->bec.m / env->consts.in.hbar;

    env->dt = 1e-4; // * (1. / sqrt(env->consts.in.G * 1 /*rho*/ ));
    env->opts.tmax = 100000;
}

//==============================================================================
//                                 write_state
//==============================================================================
void write_state(struct env *env)
{
    size_t idx=0;
    int i,j,k;

    i=0;
    if (env->space.Nz > 1) fprintf(stdout, "[");
    for (k=0; k < env->space.Nz; k++)
    {
        if (k > 0) fprintf(stdout, ",");
        if (env->space.Ny > 1) fprintf(stdout, "[");
        for (j=0; j < env->space.Ny; j++)
        {
            if (j > 0) fprintf(stdout, ",");
            fprintf(stdout, "[");
            for (i=0; i < env->space.Nx; i++)
            {
                if (i > 0) fprintf(stdout, ", ");
                //fprintf(stdout, "%.3g+%.3gj", s->phi[idx][0], s->phi[idx][1]);
                fprintf(stdout, "%.8g+%.8gj", env->state.psi[idx][0], env->state.psi[idx][1]);
                //fprintf(stdout, "%.3g+%.3gj", s->rho[idx][0], s->rho[idx][1]);
                idx++;
            }
            fprintf(stdout, "]");
        }
        if (env->space.Ny > 1) fprintf(stdout, "]");
    }
    if (env->space.Nz > 1) fprintf(stdout, "]");
    fprintf(stdout, "\n");
    fflush(stdout);
}

//==============================================================================
//                           phi_harmonic_oscillator
//==============================================================================
void phi_harmonic_oscillator(struct state *s, struct space *sp,
                             double x, double y, double z, 
                             double rx, double ry, double rz)
{
    int32_t i,j,k;
    size_t idx=0;

    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->dx + sp->xmin;
        double y0 = j*sp->dx + sp->ymin;
        double z0 = k*sp->dx + sp->zmin;

        //s->phi[idx][0] = 1000*pow(x0-x, 2);

#if 1
        if (x-rx <= x0 && x0 <= x+rx) {
        if (y-ry <= y0 && y0 <= y+ry) {
        if (z-rz <= z0 && z0 <= z+rz)
        {
            s->phi[idx][0] += 0;
            s->phi[idx][1] += 0;
        }}}
        else
        {
            s->phi[idx][0] += 1e12;
            s->phi[idx][1] += 0;
        }
#endif

        idx++;
    }
}

inline void kxyz(int i, struct space *s, double *x, double *y, double *z)
{
    int iz = i / (s->Nx * s->Ny); i -= iz * (s->Nx * s->Ny);
    int iy = i / s->Nx;           i -= iy * s->Nx;
    int ix = i;
    *x = (s->dk * ix);
    *y = (s->dk * iy);
    *z = (s->dk * iz);
}

//==============================================================================
//                                    drift
//==============================================================================
void drift(const double dt, struct env * const env, const struct plans *plans)
{
    size_t idx;
    int32_t i,j,k;

#if 0
    for (idx=0; idx < env->state.N; idx++)
    {
        //assert(fabs(env->state.psi[idx][1]) < 1e-12);
        assert(fabs(env->state.psi[idx][0] - sqrt(env->state.rho[idx][0])) < 1e-12);
    }
#endif

    fftw_execute(plans->psi_f);
    #pragma omp parallel for private(j,i,idx)
    for (k=0; k < env->space.Nz; k++) { idx = k * env->space.Nx * env->space.Ny;
    for (j=0; j < env->space.Ny; j++)
    for (i=0; i < env->space.Nx; i++)
    {
        double k2;

        const int32_t io = i - env->space.Nx*(i > (env->space.Nx>>1));
        const int32_t jo = j - env->space.Ny*(j > (env->space.Ny>>1));
        const int32_t ko = k - env->space.Nz*(k > (env->space.Nz>>1));

        //fprintf(stderr, "%ld %i %i %i\n", idx, io, jo, ko);

        k2 = pow(io*env->space.dk,2)
           + pow(jo*env->space.dk,2)
           + pow(ko*env->space.dk,2);

        //k2 /= 2;

        //fprintf(stderr, "%ld %f\n", idx, env->state.psi[idx][1]);

        const long double c = cosl(-env->drift_exp * k2 * dt);
        const long double s = sinl(-env->drift_exp * k2 * dt);
        const long double p0 = env->state.psi[idx][0];  /// sqrt(env->state.N);
        const long double p1 = env->state.psi[idx][1]; // / sqrt(env->state.N);

        env->state.psi[idx][0] = c*p0 - s*p1;
        env->state.psi[idx][1] = c*p1 + s*p0;

        //fprintf(stderr, "%.20Lf %.20Lf\n", p0*p0+p1*p1, powl(env->state.psi[idx][0],2) + powl(env->state.psi[idx][1],2));
        //assert(fdim(p0*p0+p1*p1, pow(env->state.psi[idx][0],2) + pow(env->state.psi[idx][1],2)) < 1e-12);

        idx++;
    }}

    double rho_avg = 0;
    env->rho_max = 0;
    env->rho_min = FLT_MAX;

    fftw_execute(plans->psi_b);
    #pragma omp parallel
    {
        double rho_max = 0;
        double rho_min = FLT_MAX;

        #pragma omp for reduction(+:rho_avg)
        for (idx=0; idx < env->state.N; idx++)
        {
            env->state.psi[idx][0] /= (env->state.N);
            env->state.psi[idx][1] /= (env->state.N);

            //double mag2 = pow(env->state.psi[idx][0],2) + pow(env->state.psi[idx][1],2);
            //fprintf(stderr, "%ld %.20f %.20f\n", idx, mag2, env->state.rho[idx][0]);
            //assert(fabs(env->state.psi[idx][1]) < 1e-12);
            
            //assert(fdim(mag2, env->state.rho[idx][0]) < 1e-8);

            double q = cmag2(env->state.psi[idx]);

            rho_max = MAX(rho_max, q);
            rho_min = MIN(rho_min, q);
            rho_avg += q;
        }

        #pragma omp critical
        {
            env->rho_max = MAX(env->rho_max, rho_max);
            env->rho_min = MIN(env->rho_min, rho_min);
        }
    }

    env->rho_avg = rho_avg / env->state.N;

    idx=0;
    //fprintf(stderr, "%f %f\n", env->state.psi[0][0], env->state.rho[0][0]);
}

//==============================================================================
//                                     kick
//==============================================================================
void kick(const double dt, const struct env *env)
{
    size_t idx;

    #pragma omp parallel for 
    for (idx=0; idx < env->state.N; idx++)
    {
        const long double p0  = env->state.psi[idx][0];
        const long double p1  = env->state.psi[idx][1];
        const long double phi = env->state.phi[idx][0];
        const long double c = cosl(phi * dt * env->kick_exp);
        const long double s = sinl(phi * dt * env->kick_exp);

        env->state.psi[idx][0] = c*p0 - s*p1;
        env->state.psi[idx][1] = c*p1 + s*p0;

        //fprintf(stderr, "p0/p1  %f %f -> %f %f / %f\n", p0, p1, state.psi[i][0], state.psi[i][1], state.phi[i][0]*dt);
    }
}

//==============================================================================
//                                   gravity
//==============================================================================
void gravity_periodic(struct env * const env, const struct plans *plans)
{
    size_t idx, idx_center;
    int32_t i,j,k;

    double prob_tot = 0;
    #pragma omp parallel for reduction(+:prob_tot)
    for (idx=0; idx < env->state.N; idx++)
    {
        env->state.rho[idx][0] = cmag2(env->state.psi[idx]);
        env->state.rho[idx][1] = 0;

        env->state.phi[idx][0] = 0;
        env->state.phi[idx][1] = 0;
        prob_tot += env->state.rho[idx][0];
    }

#if 0
    if (step % 100 == 0)
    {
        capture_image(&space, &state, &fb);
        write_image("frames/becon.%05i.phi.png", step, &fb);
    }
#endif

#if 1
    idx_center = env->state.N;
    idx_center >>= env->space.Nx > 1;
    idx_center >>= env->space.Ny > 1;
    idx_center >>= env->space.Nz > 1;

    Log("Total Probability is %f\n", prob_tot);
    fftw_execute(plans->rho_f);
    #pragma omp parallel for private(j,i,idx)
    for (k=0; k < env->space.Nz; k++) { idx = k * env->space.Nx * env->space.Ny;
    for (j=0; j < env->space.Ny; j++)
    for (i=0; i < env->space.Nx; i++)
    {
        double k2;

        const int32_t io = i - env->space.Nx*(i > (env->space.Nx>>1));
        const int32_t jo = j - env->space.Ny*(j > (env->space.Ny>>1));
        const int32_t ko = k - env->space.Nz*(k > (env->space.Nz>>1));

        k2 = pow(io*env->space.dk,2)
           + pow(jo*env->space.dk,2)
           + pow(ko*env->space.dk,2);

        if (k2 != 0)
        {
            env->state.phi[idx][0] = env->state.rho[idx][0] / k2;
            env->state.phi[idx][1] = env->state.rho[idx][1] / k2;
            //fprintf(stderr, "phi %g\n", env->state.phi[idx][0]);
        }
        else
        {
            /* Treating k2==0 by setting phi=0 subtracts off the mean density */
            env->state.phi[idx][0] = 0;
            env->state.phi[idx][1] = 0;
        }


        idx++;
    }}
    fftw_execute(plans->phi_b);

    env->phi_max = 0;
    #pragma omp parallel 
    {
        double phi_max = 0;
        #pragma omp for 
        for (idx=0; idx < env->state.N; idx++)
        {
            env->state.phi[idx][0] /= (env->state.N);
            env->state.phi[idx][1] /= (env->state.N);
            phi_max = MAX(phi_max, env->state.phi[idx][0]);
        }

        #pragma omp critical
        {
            env->phi_max = MAX(env->phi_max, phi_max);
        }
    }
#endif
}

void gravity_vacuum(struct env * const env, const struct plans *plans)
{
    size_t idx, vidx;
    int32_t i,j,k;

    double prob_tot = 0;
    memset(env->state.rho, 0, sizeof(*env->state.rho) * env->state.vN);

    #pragma omp parallel for reduction(+:prob_tot) private(j,i,idx, vidx)
    for (k=0; k < env->space.Nz; k++) {  idx = k * env->space.Nx * env->space.Ny;
    for (j=0; j < env->space.Ny; j++) { vidx = k * env->space.vNx * env->space.vNy + j * env->space.vNx;
    for (i=0; i < env->space.Nx; i++)
    {
        env->state.rho[vidx][0] = cmag2(env->state.psi[idx]);
        env->state.rho[vidx][1] = 0;

        prob_tot += env->state.rho[vidx][0];

        vidx++;
        idx++;
    }}}

//  struct frame_buffer fb;
//  alloc_frame_buffer(&fb, imin(256, env->space.vNx), imin(256, env->space.vNy));
//  capture_image_log(10, 12, cmap_tipsy, env, &fb);
//  write_image(&fb, "/tmp/becon.%05i.rho.png", 10);

//  static int first_time = 1;
//  if (!first_time) exit(2);
//  first_time = 0;
    

#if 0
    if (step % 100 == 0)
    {
        capture_image(&space, &state, &fb);
        write_image("frames/becon.%05i.phi.png", step, &fb);
    }
#endif

    //Log("Total Probability is %f\n", prob_tot);
    fftw_execute(plans->rho_f);

    #pragma omp parallel for private(j,i,vidx)
    for (k=0; k < env->space.vNz; k++) { vidx = k * env->space.vNx * env->space.vNy;
    for (j=0; j < env->space.vNy; j++)
    for (i=0; i < env->space.vNx; i++)
    {
        env->state.rho[vidx][0] *= env->state.S[vidx][0];
        env->state.rho[vidx][1] *= env->state.S[vidx][1];

        vidx++;
    }}

    fftw_execute(plans->rho_b);

    env->phi_max = 0;
    {
        double phi_max = 0;
        #pragma omp for private(j,i,idx,vidx)
        for (k=0; k < env->space.Nz; k++) {  idx = k * env->space.Nx * env->space.Ny;
        for (j=0; j < env->space.Ny; j++) { vidx = k * env->space.vNx * env->space.vNy + j * env->space.vNx;
        for (i=0; i < env->space.Nx; i++)
        {
            env->state.phi[idx][0] = env->state.rho[vidx][0] / env->state.vN;
            env->state.phi[idx][1] = env->state.rho[vidx][1] / env->state.vN;
            //phi_max = MAX(phi_max, env->state.phi[idx][0]);

            idx++;
            vidx++;
        }}}

        //#pragma omp critical
        {
            //env->phi_max = MAX(env->phi_max, phi_max);
        }
    }
}

//==============================================================================
//                                     main
//==============================================================================
int main(int argc, char **argv)
{
    size_t idx;
    int32_t i,j,k;
    size_t nthreads = omp_get_num_threads();
    struct env env;
    struct plans plans;
    struct frame_buffer fb;
    struct matrix rrho;
    struct matrix rphi;

    double t;
    int step;

    void (*gravity)(struct env * const, const struct plans *) = NULL;

    env.opts.with_gravity = 1;
    env.opts.with_vacuum = 0;
    //env.opts.dt = 0.01 * 0.861522;
    //env.opts.dt = 0.01 * 0.094638;
    env.opts.dt = 1.0 * 0.013410;
    env.opts.tmax = env.opts.dt * 10000;

    env.opts.tmax = 100000;

    env.cosmo.a_start = 1;
    env.cosmo.H0 = 0.70;
    env.cosmo.omega_m = 0.279;
    env.cosmo.omega_v = 0.721;
    env.cosmo.omega_r = 0;
    env.cosmo.omega_k = 1 - (env.cosmo.omega_m + env.cosmo.omega_v + env.cosmo.omega_r);

    env.consts.si.c = 3e8;  // [m/s]
    env.consts.si.hbar = 6.58211928e-16;  // [eV s]
    env.consts.si.H0 = env.cosmo.H0 / 3.08568025e19; // [s^-1]
    env.consts.si.G  = 6.6738480e-11; // [m^3 kg^-1 s^-2]

    env.bec.m = 2.5e-22 / pow(env.consts.si.c,2);   // [eV / c^2]
    env.bec.m = 2e3 / pow(env.consts.si.c,2);   // [eV / c^2]

#if 1
    if (alloc_grafic("graficICs/128/level0", &env) != 0)
    {
        Log("Failed to read input.\n");
        exit(1);
    }
#endif

    //ic_init_test1D(&env);
    //ic_init_spherical_collapse(&env);
    //ic_init_infinite_sheet(&env);

    #pragma omp parallel for 
    for (idx=0; idx < env.state.N; idx++)
    {
        env.state.psi[idx][0] = 0;
        env.state.psi[idx][1] = 0;

        env.state.phi[idx][0] = 0;
        env.state.phi[idx][1] = 0;

        //env.state.rho[idx][0] = 0;
        //env.state.rho[idx][1] = 0;
    }

    alloc_frame_buffer(&fb, imin(256, env.space.Nx), imin(256, env.space.Ny));
    alloc_matrix(&rrho, env.space.Nmax/2+1, 3);
    alloc_matrix(&rphi, env.space.Nmax/2+1, 3);

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);


    Log("Creating plans...\n");
#define PLAN(var, dir) \
    fftw_plan_dft_3d(env.space.Nx, env.space.Ny, env.space.Nz,env.state. var,env.state. var, dir, FFTW_MEASURE)
#define vPLAN(var, dir) \
    fftw_plan_dft_3d(env.space.vNx, env.space.vNy, env.space.vNz,env.state. var,env.state. var, dir, FFTW_MEASURE)

    plans.psi_f = PLAN(psi, FFTW_FORWARD);
    plans.phi_f = PLAN(phi, FFTW_FORWARD);
    plans.psi_b = PLAN(psi, FFTW_BACKWARD);
    plans.phi_b = PLAN(phi, FFTW_BACKWARD);

    if (env.opts.with_vacuum)
    {
        plans.rho_f  = vPLAN(rho, FFTW_FORWARD);
        plans.rho_b  = vPLAN(rho, FFTW_BACKWARD);
        plans.vphi_b = vPLAN(vphi, FFTW_BACKWARD);
        plans.S_f    = vPLAN(S, FFTW_FORWARD);
        //plans.vphi_b = fftw_plan_dft_3d(env.space.vNx, env.space.vNy, env.space.vNz, env.state.vphi,env.state.vphi, FFTW_BACKWARD, FFTW_MEASURE);
        //plans.S_f    = fftw_plan_dft_3d(env.space.vNx, env.space.vNy, env.space.vNz, env.state.S,   env.state.S,    FFTW_FORWARD,  FFTW_MEASURE);
    }
    else
    {
        plans.rho_f = PLAN(rho, FFTW_FORWARD);
        plans.rho_b = PLAN(rho, FFTW_BACKWARD);
    }

#if 0
    double a = env.cosmo.a_start;
    double H2 = pow(env.cosmo.H0, 2)
             * (env.cosmo.omega_v + ((((env.cosmo.omega_r / a) + env.cosmo.omega_m) / a + env.cosmo.omega_k) / a / a));

    //XXX
    env.cosmo.rho_crit = 3*H2 / (8*M_PI);
#endif

#if 1
    if (read_grafic("graficICs/128/level0", &env) != 0)
    {
        Log("Failed to load input.\n");
        exit(1);
    }
#endif

    //ic_test1D(&env);
    //ic_spherical_collapse(&env);
    //ic_infinite_sheet(&env);

    if (env.opts.with_vacuum)
    {
        for (k=0; k < env.space.vNz; k++) { idx = k * env.space.vNx * env.space.vNy;
        for (j=0; j < env.space.vNy; j++) 
        for (i=0; i < env.space.vNx; i++)
        {
            const int32_t io = i;// - env.space.vNx*(i > (env.space.vNx>>1));
            const int32_t jo = j;// - env.space.vNy*(j > (env.space.vNy>>1));
            const int32_t ko = k;// - env.space.vNz*(k > (env.space.vNz>>1));

            env.state.S[idx][0] = 1.0;
#if 0
            / sqrt(pow(io*env.space.dx,2)

                                           + pow(jo*env.space.dx,2)
                                           + pow(ko*env.space.dx,2) 
                                           + pow(env.space.dx,2));
#endif
            env.state.S[idx][1] = 0;

            idx++;
        }}

        fftw_execute(plans.S_f);
    }

    //env.cosmo.m = 1;
    //double v = sqrt(env.cosmo.rho_crit) * pow(env.space.dx,2) * env.cosmo.m;
    //env.cosmo.h = 1e2; //300 * v;
    //env.cosmo.h_m = env.cosmo.h / env.cosmo.m;
    env.cosmo.a = env.cosmo.a_start;

    //env.eta = 1.119449183; /* XXX: only valid for 128^3 sim */


    Log("rho_crit at a=%0.4g is %0.4g\n", env.cosmo.a_start, env.cosmo.rho_crit);
    Log("bec m      = %g\n", env.bec.m);
    Log("eta        = %g\n", env.eta);
    //Log("hbar     = %g\n", env.cosmo.h);
    //Log("m        = %g\n", env.cosmo.m);
    //Log("hbar / m = %g\n", env.cosmo.h_m);
    //Log("LambdaDB = %g\n", env.cosmo.h_m / v);


    assert(env.space.Nx == 1 || (env.space.Nx & 1) == 0);
    assert(env.space.Ny == 1 || (env.space.Ny & 1) == 0);
    assert(env.space.Nz == 1 || (env.space.Nz & 1) == 0);

    Log("Running simulation...\n");

    //write_state(&env);

    //capture_image(&space, &state, &fb);
    //write_image("frames/becon.%05i.phi.png", 0, &fb);

    capture_image_log(-5, 5, cmap_tipsy, &env, &fb);
    write_image(&fb, "/tmp/becon.%05i.rho.png", 0);
    //write_matrix(&rrho, "/tmp/becon.%05i.rho.txt", 0);

    //goto shutdown;

    //write_tipsy_grid(&env, "/tmp/becon.grid.bin");


    if (env.opts.with_vacuum)
    {
        gravity = gravity_vacuum;
    }
    else
    {
        gravity = gravity_periodic;
    }


    if (env.opts.with_gravity)
    {
        gravity(&env, &plans);
    }

    //env.opts.dt = env.dt * 100;
    env.dt = env.opts.dt;
//  capture_image_log(-12, 12, cmap_tipsy, &env, &fb);
//  write_image(&fb, "/tmp/becon.%05i.rho.png", 1);
//  exit(1);

    double z = 1/env.cosmo.a - 1;
    double dz = z / 100;
    double kmax2 = pow(env.space.Nx,2) + pow(env.space.Ny,2) + pow(env.space.Nz,2);
    //for (t=0, step=0; z >= 0; z -= dz) 
    for (t=0, step=0; t < env.opts.tmax; t += env.dt)
    {
        env.cosmo.a = 1; // / (z+1);

        Log("Time %8.4g  Step %i  dt %8.4g  z %8.4f\n", t, step+1, env.dt, z);

        /* XXX: Be careful here to change drift_exp at the right time.
         * Shouldn't do this before a complete timestep is finished */
        env.drift_exp = 1. / (2 * env.eta * pow(env.cosmo.a, 2));
        env.kick_exp  = 3 * env.cosmo.omega_m * env.eta / (2 * env.cosmo.a);

        //env.dt = fabs(M_PI * env.cosmo.a / (6 * env.cosmo.omega_m * env.eta * env.phi_max));
        //

        //fprintf(stderr, "%f %f\n", env.dt, env.kick_exp * M_PI_4);
        //fprintf(stderr, "%f %f\n", env.dt, env.drift_exp * M_PI_4);

        /* Check time step upper bound */
        assert(env.kmax2 * env.drift_exp * env.dt/2 <= M_PI_4);

        /* Check time step lower bound */
        //assert(env.kick_exp * env.dt <= M_PI_4);

        drift(env.dt/2, &env, &plans);

        if (env.opts.with_gravity)
        {
            gravity(&env, &plans);
        }

        //t += env.dt;

        if (t >= env.opts.dt * (step+1))
        {
            //write_state(&env);

            step++;
            t = step * env.opts.dt;

#if 1
            Log("rho_min is %.8g\n", env.rho_min);
            Log("rho_max is %.8g\n", env.rho_max);
            Log("rho_avg is %.8g\n", env.rho_avg);
            //capture_image_log(log10(env.rho_max)*0.8, 1.1*log10(env.rho_max), cmap_tipsy, &env, &fb);
            capture_image_log(-1, 1, cmap_tipsy, &env, &fb);
            //capture_image_log(1, 10*env.rho_max, cmap_tipsy, &env, &fb);
            //capture_image(&space, &state, &fb);
            write_image(&fb, "/tmp/becon.%05i.rho.png", step);
            //write_ascii_rho(&env, "/tmp/becon.%05i.rho", step);

//          radial_avg_density(&env, &rrho, 0,0,0);
//          radial_avg_potential(&env, &rphi, 0,0,0);
//          write_matrix(&rrho, "/tmp/becon.%05i.rho.txt", step);
//          write_matrix(&rphi, "/tmp/becon.%05i.phi.txt", step);
#endif
        }

        kick(env.dt, &env);


        //phi_harmonic_oscillator(&state, &space, 0,0,0, space.xmax/2, 0,0);

#if 0
        fftw_execute(plans.rho_b);
        for (idx=0; idx < state.N; idx++)
        {
            state.rho[idx][0] /= state.N;
            state.rho[idx][1] /= state.N;
        }
#endif

    }

shutdown:

    Log("Simulation complete.\n");


    fftw_destroy_plan(plans.psi_f);
    fftw_destroy_plan(plans.phi_f);
    fftw_destroy_plan(plans.rho_f);

    fftw_destroy_plan(plans.psi_b);
    fftw_destroy_plan(plans.phi_b);
    fftw_destroy_plan(plans.rho_b);
    fftw_cleanup_threads();

    fftw_free(env.state.psi);
    fftw_free(env.state.phi);
    fftw_free(env.state.rho);

    if (env.opts.with_vacuum)
    {
        fftw_free(env.state.S);
        fftw_free(env.state.vphi);
    }

    free_frame_buffer(&fb);

    return EXIT_SUCCESS;
}

