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

    Log("    %ld grid cells.\n", N);

    state->psi = (fftw_complex*) fftw_malloc(sizeof(*state->psi) * N);
    state->phi = (fftw_complex*) fftw_malloc(sizeof(*state->phi) * N);
    state->rho = (fftw_complex*) fftw_malloc(sizeof(*state->rho) * N);
}

void ic_test1D(struct state *s, struct space *sp)
{
    int32_t x;

    Log("Creating test1D initial conditions...\n");

    for (x=0; x < sp->Nx; x++)
    {
        s->psi[x][0] = exp(-140 * pow(sp->xmin + x*sp->drx, 2));

        //s->psi[x][0] = cos(2*M_PI * x/sp->Nx);
        //s->psi[x][0] = 1e-5 * cos(2*M_PI * x/L);
        s->psi[x][1] = 0; //sin(2*M_PI * x/L);

        s->phi[x][0] = 0;
        s->phi[x][1] = 0;

        s->rho[x][0] = 0;
        s->rho[x][1] = 0;
    }
}

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
                i++;
            }
            fprintf(stdout, "]");
        }
        if (sp->Ny > 1) fprintf(stdout, "]");
    }
    if (sp->Nz > 1) fprintf(stdout, "]");
    fprintf(stdout, "\n");
    fflush(stdout);
}

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

int main(int argc, char **argv)
{
    size_t idx;
    size_t nthreads = 2; //omp_get_num_threads();
    struct state state;
    struct space space;
    struct plans plans;
    struct options opts;
    struct frame_buffer fb;

    int32_t i,j,k;

    //double dt = 1e-8;
    //double tmax = dt * 0; //1.0001;; //1.0;
    //double tmax = dt * 10000; //1.0001;; //1.0;
    double t;
    int step;

    opts.with_gravity = 1;
    opts.dt = 1e-5;
    opts.tmax = opts.dt * 1000;

    if (alloc_grafic("graficICs/level0", &state, &space) != 0)
    {
        Log("Failed to read input.\n");
        exit(1);
    }

    alloc_frame_buffer(&fb, imin(256, space.Nx), imin(256, space.Ny));

    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);

    //ic_init_test1D(&state, &space);

    Log("Creating plans...\n");
#define PLAN(sp, st, var, dir) \
    fftw_plan_dft_3d(sp.Nx, sp.Ny, sp.Nz, st. var, st. var, dir, FFTW_ESTIMATE)

    plans.psi_f = PLAN(space, state, psi, FFTW_FORWARD);
    plans.phi_f = PLAN(space, state, phi, FFTW_FORWARD);
    plans.rho_f = PLAN(space, state, rho, FFTW_FORWARD);
    plans.psi_b = PLAN(space, state, psi, FFTW_BACKWARD);
    plans.phi_b = PLAN(space, state, phi, FFTW_BACKWARD);
    plans.rho_b = PLAN(space, state, rho, FFTW_BACKWARD);

    if (read_grafic("graficICs/level0", &state, &space) != 0)
    {
        Log("Failed to load input.\n");
        exit(1);
    }

    //ic_test1D(&state, &space);

    assert((space.Nx & 1) == 0);
    assert((space.Ny & 1) == 0);
    assert((space.Nz & 1) == 0);

    Log("Running simulation...\n");

    //write_state(&state, &space);

    for (t=0, step=0; t < opts.tmax; t += opts.dt, step++)
    {
        Log("Time %8.4g  Step %i\n", t, step);

        //----------------------------------------------------------------------
        // Drift
        //----------------------------------------------------------------------
        fftw_execute(plans.psi_f);
        size_t idx = 0;
        for (k=0; k < space.Nz; k++)
        for (j=0; j < space.Ny; j++)
        for (i=0; i < space.Nx; i++)
        {
            double k2, x,y,z;

            //kxyz(i, &space, &x,&y,&z);
            //k2 = x*x + y*y + z*z;

            k2 = pow(((i + (space.Nx>>1))%space.Nx)*space.dkx,2)
               + pow(((j + (space.Ny>>1))%space.Ny)*space.dky,2)
               + pow(((k + (space.Nz>>1))%space.Nz)*space.dkz,2);

            k2 /= 2;

            double c = cos(k2 * opts.dt);
            double s = sin(k2 * opts.dt);
            double p0 = state.psi[idx][0];
            double p1 = state.psi[idx][1];

            state.psi[idx][0] = c*p0 - s*p1;
            state.psi[idx][1] = c*p1 + s*p0;

            idx++;
        }

        fftw_execute(plans.psi_b);
        for (idx=0; idx < state.N; idx++)
        {
            state.psi[idx][0] /= state.N;
            state.psi[idx][1] /= state.N;
        }

        if (opts.with_gravity)
        {
            //----------------------------------------------------------------------
            // Gravity
            //----------------------------------------------------------------------
            double prob_tot = 0;
            for (idx=0; idx < state.N; idx++)
            {
                state.rho[idx][0] = pow(state.psi[idx][0],2) + pow(state.psi[idx][1],2);
                state.rho[idx][1] = 0;

                state.phi[idx][0] = 0;
                state.phi[idx][1] = 0;
                prob_tot += state.rho[idx][0];
            }

            //Log("Total Probability is %f\n", prob_tot);
            fftw_execute(plans.rho_f);
            idx = 0;
            for (k=0; k < space.Nz; k++)
            for (j=0; j < space.Ny; j++)
            for (i=0; i < space.Nx; i++)
            {
                double k2, x,y,z;

                //kxyz(i, &space, &x,&y,&z);
                //k2 = x*x + y*y + z*z;

                k2 = pow(((i + (space.Nx>>1))%space.Nx)*space.dkx,2)
                   + pow(((j + (space.Ny>>1))%space.Ny)*space.dky,2)
                   + pow(((k + (space.Nz>>1))%space.Nz)*space.dkz,2);

                if (fabs(k2) > 1e-12)
                {
                    state.phi[idx][0] = state.rho[idx][0] / k2;
                    state.phi[idx][1] = 0;
                }

                idx++;
            }
            fftw_execute(plans.phi_b);
            for (idx=0; idx < state.N; idx++)
            {
                state.phi[idx][0] /= state.N;
                state.phi[idx][1] /= state.N;
            }
        }

        phi_harmonic_oscillator(&state, &space, 0,0,0, space.xmax/2, 0,0);

        //----------------------------------------------------------------------
        // Kick
        //----------------------------------------------------------------------
        for (idx=0; idx < state.N; idx++)
        {
            double p0 = state.psi[idx][0];
            double p1 = state.psi[idx][1];
            double c = cos(state.phi[idx][0] * opts.dt/2.);
            double s = sin(state.phi[idx][0] * opts.dt/2.);

            state.psi[idx][0] = c*p0 - s*p1;
            state.psi[idx][1] = c*p1 + s*p0;

            //fprintf(stderr, "p0/p1  %f %f -> %f %f / %f\n", p0, p1, state.psi[i][0], state.psi[i][1], state.phi[i][0]*dt);
        }

#if 1
        fftw_execute(plans.rho_b);
        for (idx=0; idx < state.N; idx++)
        {
            state.rho[idx][0] /= state.N;
            state.rho[idx][1] /= state.N;
        }
#endif

        if (step % 10 == 0)
        {
            //write_state(&state, &space);
            capture_image(&space, &state, &fb);
            write_image("frames/becon.%05i.jpg", step, &fb);
        }
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

