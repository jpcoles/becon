#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fftw3.h>

#include "becon.h"
#include "log.h"
#include "io_grafic.h"
#include "io_image.h"
#include "frame_buffer.h"
#include "cmap.h"

inline int imin(int a, int b)
{
    if (a < b) return a;
    return b;
}

inline int imax(int a, int b)
{
    if (a > b) return a;
    return b;
}


//==============================================================================
//                                ic_init_test1D
//==============================================================================
void ic_init_test1D(struct state *state, struct space *space)
{
    Log("Preparing test1D initial conditions...\n");

    space->xmin = -1.5;
    space->xmax = 1.5;
    space->ymin = 0.0;
    space->ymax = 0.0;
    space->zmin = 0.0;
    space->zmax = 0.0;

    space->Nx = 1 << 12;
    space->Ny = 1;
    space->Nz = 1;

    space->drx = (space->xmax-space->xmin) / space->Nx;
    space->dry = (space->ymax-space->ymin) / space->Ny;
    space->drz = (space->zmax-space->zmin) / space->Nz;

    space->dkx = fabs(space->drx) != 0 ? 2*M_PI/(space->drx * space->Nx) : 0;
    space->dky = fabs(space->dry) != 0 ? 2*M_PI/(space->dry * space->Ny) : 0;
    space->dkz = fabs(space->drz) != 0 ? 2*M_PI/(space->drz * space->Nz) : 0;

    //Log("dt: %g\n", dt);
    Log("    dx: %g %g %g\n", space->drx, space->dry, space->drz);
    Log("    dk: %g %g %g\n", space->dkx, space->dky, space->dkz);

    size_t N = space->Nx * space->Ny * space->Nz;

    state->N = N;

    Log("    %ld grid cells.\n", N);

    state->psi = (fftw_complex*) fftw_malloc(sizeof(*state->psi) * N);
    state->phi = (fftw_complex*) fftw_malloc(sizeof(*state->phi) * N);
    state->rho = (fftw_complex*) fftw_malloc(sizeof(*state->rho) * N);
}

//==============================================================================
//                                  ic_test1D
//==============================================================================
void ic_test1D(struct state *s, struct space *sp)
{
    int32_t x;

    Log("Creating test1D initial conditions...\n");

    for (x=0; x < sp->Nx; x++)
    {
        s->psi[x][0] = exp(-1 * pow(sp->xmin + x*sp->drx, 2));
        //fprintf(stderr, "%f\n", s->psi[x][0]);

        //s->psi[x][0] = cos(2*M_PI * x/sp->Nx);
        //s->psi[x][0] = 1e-5 * cos(2*M_PI * x/L);
        s->psi[x][1] = 0; //sin(2*M_PI * x/L);

        s->phi[x][0] = 0;
        s->phi[x][1] = 0;

        s->rho[x][0] = 0;
        s->rho[x][1] = 0;
    }
}

//==============================================================================
//                                 write_state
//==============================================================================
void write_state(struct state *s, struct space *sp)
{
    size_t idx=0;
    int i,j,k;

    i=0;
    if (sp->Nz > 1) fprintf(stdout, "[");
    for (k=0; k < sp->Nz; k++)
    {
        if (k > 0) fprintf(stdout, ",");
        if (sp->Ny > 1) fprintf(stdout, "[");
        for (j=0; j < sp->Ny; j++)
        {
            if (j > 0) fprintf(stdout, ",");
            fprintf(stdout, "[");
            for (i=0; i < sp->Nx; i++)
            {
                if (i > 0) fprintf(stdout, ", ");
                //fprintf(stdout, "%.3g+%.3gj", s->phi[idx][0], s->phi[idx][1]);
                fprintf(stdout, "%.3g+%.3gj", s->psi[idx][0], s->psi[idx][1]);
                //fprintf(stdout, "%.3g+%.3gj", s->rho[idx][0], s->rho[idx][1]);
                idx++;
            }
            fprintf(stdout, "]");
        }
        if (sp->Ny > 1) fprintf(stdout, "]");
    }
    if (sp->Nz > 1) fprintf(stdout, "]");
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
        double x0 = i*sp->drx + sp->xmin;
        double y0 = j*sp->dry + sp->ymin;
        double z0 = k*sp->drz + sp->zmin;

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
    *x = (s->dkx * ix);
    *y = (s->dky * iy);
    *z = (s->dkz * iz);
}

//==============================================================================
//                                    drift
//==============================================================================
void drift(const double dt, const struct space *space, const struct state *state, const struct plans *plans)
{
    size_t idx;
    int32_t i,j,k;

#if 0
    for (idx=0; idx < state->N; idx++)
    {
        //assert(fabs(state->psi[idx][1]) < 1e-12);
        assert(fabs(state->psi[idx][0] - sqrt(state->rho[idx][0])) < 1e-12);
    }
#endif

    fftw_execute(plans->psi_f);
    #pragma omp parallel for private(j,i,idx)
    for (k=0; k < space->Nz; k++) { idx = k * space->Nx * space->Ny;
    for (j=0; j < space->Ny; j++)
    for (i=0; i < space->Nx; i++)
    {
        double k2;

        const int32_t io = i - space->Nx*(i > (space->Nx>>1));
        const int32_t jo = j - space->Ny*(j > (space->Ny>>1));
        const int32_t ko = k - space->Nz*(k > (space->Nz>>1));

        //fprintf(stderr, "%ld %i %i %i\n", idx, io, jo, ko);

        k2 = pow(io*space->dkx,2)
           + pow(jo*space->dky,2)
           + pow(ko*space->dkz,2);

        k2 /= 2;

        //fprintf(stderr, "%ld %f\n", idx, state->psi[idx][1]);

        const long double c = cosl(-k2 * dt);
        const long double s = sinl(-k2 * dt);
        const long double p0 = state->psi[idx][0];// / sqrt(state->N);
        const long double p1 = state->psi[idx][1];// / sqrt(state->N);

        state->psi[idx][0] = c*p0 - s*p1;
        state->psi[idx][1] = c*p1 + s*p0;

        //fprintf(stderr, "%.20Lf %.20Lf\n", p0*p0+p1*p1, powl(state->psi[idx][0],2) + powl(state->psi[idx][1],2));
        //assert(fdim(p0*p0+p1*p1, pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2)) < 1e-12);

        idx++;
    }}

    fftw_execute(plans->psi_b);
    #pragma omp parallel for 
    for (idx=0; idx < state->N; idx++)
    {
        state->psi[idx][0] /= (state->N);
        state->psi[idx][1] /= (state->N);

        //double mag2 = pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2);
        //fprintf(stderr, "%ld %.20f %.20f\n", idx, mag2, state->rho[idx][0]);
        //assert(fabs(state->psi[idx][1]) < 1e-12);
        
        //assert(fdim(mag2, state->rho[idx][0]) < 1e-8);
    }

    idx=0;
    //fprintf(stderr, "%f %f\n", state->psi[0][0], state->rho[0][0]);
}

//==============================================================================
//                                     kick
//==============================================================================
void kick(const double dt, const struct space *space, const struct state *state, const struct plans *plans)
{
    size_t idx;

    #pragma omp parallel for 
    for (idx=0; idx < state->N; idx++)
    {
        double p0 = state->psi[idx][0];
        double p1 = state->psi[idx][1];
        double c = cos(state->phi[idx][0] * dt);
        double s = sin(state->phi[idx][0] * dt);

        state->psi[idx][0] = c*p0 - s*p1;
        state->psi[idx][1] = c*p1 + s*p0;

        //fprintf(stderr, "p0/p1  %f %f -> %f %f / %f\n", p0, p1, state.psi[i][0], state.psi[i][1], state.phi[i][0]*dt);
    }
}

//==============================================================================
//                                   gravity
//==============================================================================
void gravity(const struct space *space, const struct state *state, const struct plans *plans)
{
    size_t idx;
    int32_t i,j,k;

    double prob_tot = 0;
    #pragma omp parallel for reduction(+:prob_tot)
    for (idx=0; idx < state->N; idx++)
    {
        state->rho[idx][0] = pow(state->psi[idx][0],2) + pow(state->psi[idx][1],2);
        state->rho[idx][1] = 0;

        state->phi[idx][0] = 0;
        state->phi[idx][1] = 0;
        prob_tot += state->rho[idx][0];
    }

#if 0
    if (step % 100 == 0)
    {
        capture_image(&space, &state, &fb);
        write_image("frames/becon.%05i.phi.jpg", step, &fb);
    }
#endif

#if 1
    //Log("Total Probability is %f\n", prob_tot);
    fftw_execute(plans->rho_f);
    #pragma omp parallel for private(j,i,idx)
    for (k=0; k < space->Nz; k++) { idx = k * space->Nx * space->Ny;
    for (j=0; j < space->Ny; j++)
    for (i=0; i < space->Nx; i++)
    {
        double k2;

        const int32_t io = i <= (space->Nx>>1) ? i : i-(space->Nx);
        const int32_t jo = j <= (space->Ny>>1) ? j : j-(space->Ny);
        const int32_t ko = k <= (space->Nz>>1) ? k : k-(space->Nz);

        k2 = pow(io*space->dkx,2)
           + pow(jo*space->dky,2)
           + pow(ko*space->dkz,2);

        if (fabs(k2) > 1e-12)
        {
            state->phi[idx][0] = state->rho[idx][0] / k2; // / (state->N);
            //fprintf(stderr, "phi %g\n", state->phi[idx][0]);
            state->phi[idx][1] = 0;
        }
        else
        {
            state->phi[idx][0] = 0;
            state->phi[idx][1] = 0;
        }


        idx++;
    }}
    fftw_execute(plans->phi_b);
    #pragma omp parallel for 
    for (idx=0; idx < state->N; idx++)
    {
        state->phi[idx][0] /= (state->N);
        state->phi[idx][1] /= (state->N);
    }
#endif
}

//==============================================================================
//                                     main
//==============================================================================
int main(int argc, char **argv)
{
    size_t nthreads = 6; //omp_get_num_threads();
    struct state state;
    struct space space;
    struct plans plans;
    struct options opts;
    struct frame_buffer fb;

    double t;
    int step;

    opts.with_gravity = 1;
    opts.dt = 0.01 * 0.861522;
    opts.tmax = opts.dt * 10000;

#if 1
    if (alloc_grafic("graficICs/256/level0", &state, &space) != 0)
    {
        Log("Failed to read input.\n");
        exit(1);
    }
#endif

    //ic_init_test1D(&state, &space);

    alloc_frame_buffer(&fb, imin(256, space.Nx), imin(256, space.Ny));

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);


    Log("Creating plans...\n");
#define PLAN(sp, st, var, dir) \
    fftw_plan_dft_3d(sp.Nx, sp.Ny, sp.Nz, st. var, st. var, dir, FFTW_MEASURE)

    plans.psi_f = PLAN(space, state, psi, FFTW_FORWARD);
    plans.phi_f = PLAN(space, state, phi, FFTW_FORWARD);
    plans.rho_f = PLAN(space, state, rho, FFTW_FORWARD);
    plans.psi_b = PLAN(space, state, psi, FFTW_BACKWARD);
    plans.phi_b = PLAN(space, state, phi, FFTW_BACKWARD);
    plans.rho_b = PLAN(space, state, rho, FFTW_BACKWARD);

#if 1
    if (read_grafic("graficICs/256/level0", &state, &space) != 0)
    {
        Log("Failed to load input.\n");
        exit(1);
    }
#endif

    //ic_test1D(&state, &space);

    assert(space.Nx == 1 || (space.Nx & 1) == 0);
    assert(space.Ny == 1 || (space.Ny & 1) == 0);
    assert(space.Nz == 1 || (space.Nz & 1) == 0);

    Log("Running simulation...\n");

    //write_state(&state, &space);

    //capture_image(&space, &state, &fb);
    //write_image("frames/becon.%05i.phi.jpg", 0, &fb);

    //capture_image_log(1, 100, cmap_tipsy, &space, &state, &fb);
    //write_image("frames/becon.%05i.rho.jpg", 0, &fb);

    for (t=0, step=1; t < opts.tmax; t = step++ * opts.dt)
    {
        Log("Time %8.4g  Step %i\n", t, step);

        drift(opts.dt/2, &space, &state, &plans);

#if 1
        if (step % 10 == 0)
        {
            capture_image_log(1, 100, cmap_tipsy, &space, &state, &fb);
            //capture_image(&space, &state, &fb);
            write_image("frames/becon.%05i.rho.jpg", step, &fb);
        }
#endif

        if (opts.with_gravity)
        {
            gravity(&space, &state, &plans);
        }

        kick(opts.dt, &space, &state, &plans);


        //phi_harmonic_oscillator(&state, &space, 0,0,0, space.xmax/2, 0,0);

#if 0
        fftw_execute(plans.rho_b);
        for (idx=0; idx < state.N; idx++)
        {
            state.rho[idx][0] /= state.N;
            state.rho[idx][1] /= state.N;
        }
#endif

        //write_state(&state, &space);
    }

    Log("Simulation complete.\n");


    fftw_destroy_plan(plans.psi_f);
    fftw_destroy_plan(plans.phi_f);
    fftw_destroy_plan(plans.rho_f);

    fftw_destroy_plan(plans.psi_b);
    fftw_destroy_plan(plans.phi_b);
    fftw_destroy_plan(plans.rho_b);
    fftw_cleanup_threads();

    fftw_free(state.psi);
    fftw_free(state.phi);
    fftw_free(state.rho);

    free_frame_buffer(&fb);

    return EXIT_SUCCESS;
}

