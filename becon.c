#include <assert.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fftw3.h>

#include <omp.h>

#include "becon.h"

#define XYZ(i, x,y,z, dx,dy,dz) do { \
    int _i = i;\
    int iz = _i / (space.Nx*space.Ny); \
    _i -= iz*(space.Nx * space.Ny); \
    int iy = _i / space.Nx; \
    _i -= iy*space.Nx; \
    int ix = _i; \
    x = (dx * ix); \
    y = (dy * iy); \
    z = (dz * iz); \
} while (0)

#define XYZ0(i, x,y,z, dx,dy,dz) do { \
    int _i = i;\
    int iz = _i / (rank_size[0]*rank_size[1]); \
    _i -= iz*(rank_size[0]*rank_size[1]); \
    int iy = _i / rank_size[0]; \
    _i -= iy*rank_size[0]; \
    int ix = _i; \
    x = (dx * ix); \
    y = (dy * iy); \
    z = (dz * iz); \
} while (0)

struct state
{
    fftw_complex *psi;
    fftw_complex *phi;
    fftw_complex *rho;
};

struct space
{
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    int Nx, Ny, Nz;
    double drx, dry, drz;
    double dkx, dky, dkz;
};

struct plans
{
    fftw_plan psi_f;
    fftw_plan phi_f;
    fftw_plan rho_f;

    fftw_plan psi_b;
    fftw_plan phi_b;
    fftw_plan rho_b;
};


void Log(const char *fmt, ...)
{
    char s[100];
    va_list ap;
    time_t t;
    struct tm *tmp;
    t = time(NULL);
    tmp = localtime(&t);

    strftime(s, sizeof(s), "%F %T", tmp);

    va_start(ap, fmt);
    fprintf(stderr, "%15s  |  ", s);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
}

void ic_test1D(struct state *s, struct space *sp)
{
    int x;

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
    int i,x,y,z;

    i=0;
    if (sp->Nz > 1) fprintf(stdout, "[");
    for (z=0; z < sp->Nz; z++)
    {
        if (z > 0) fprintf(stdout, ",");
        if (sp->Ny > 1) fprintf(stdout, "[");
        for (y=0; y < sp->Ny; y++)
        {
            if (y > 0) fprintf(stdout, ",");
            fprintf(stdout, "[");
            for (x=0; x < sp->Nx; x++)
            {
                if (x > 0) fprintf(stdout, ", ");
                //fprintf(stdout, "%.3g+%.3gj", s->phi[i][0], s->phi[i][1]);
                fprintf(stdout, "%.3g+%.3gj", s->psi[i][0], s->psi[i][1]);
                //fprintf(stdout, "%.3g+%.3gj", s->rho[i][0], s->rho[i][1]);
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
    int i,j,k;

    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->drx + sp->xmin;
        double y0 = j*sp->dry + sp->ymin;
        double z0 = k*sp->drz + sp->zmin;

        int I = k * (sp->Nx*sp->Ny) + j*sp->Nx + i;
        //s->phi[I][0] = 1000*pow(x0-x, 2);

#if 1
        if (x-rx <= x0 && x0 <= x+rx) {
        if (y-ry <= y0 && y0 <= y+ry) {
        if (z-rz <= z0 && z0 <= z+rz)
        {
            int I = k * (sp->Nx*sp->Ny) + j*sp->Nx + i;
            s->phi[I][0] += 0;
            s->phi[I][1] += 0;
        }}}
        else
        {
            int I = k * (sp->Nx*sp->Ny) + j*sp->Nx + i;
            s->phi[I][0] += 1e12;
            s->phi[I][1] += 0;
        }
#endif
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
    int i;
    size_t nthreads = 1; //omp_get_num_threads();
    struct state state;
    struct space space;
    struct plans plans;

    double dt = 1e-8;
    double tmax = dt * 10000; //1.0001;; //1.0;
    double t;
    int step;

    space.xmin = -1.5;
    space.xmax = 1.5;
    space.ymin = 0.0;
    space.ymax = 0.0;
    space.zmin = 0.0;
    space.zmax = 0.0;

    space.Nx = 1 << 12;
    space.Ny = 1;
    space.Nz = 1;

    space.drx = (space.xmax-space.xmin) / space.Nx;
    space.dry = (space.ymax-space.ymin) / space.Ny;
    space.drz = (space.zmax-space.zmin) / space.Nz;

    space.dkx = fabs(space.drx) != 0 ? 2*M_PI/(space.drx * space.Nx) : 0;
    space.dky = fabs(space.dry) != 0 ? 2*M_PI/(space.dry * space.Ny) : 0;
    space.dkz = fabs(space.drz) != 0 ? 2*M_PI/(space.drz * space.Nz) : 0;

    Log("dx: %f %f %f\n", space.drx, space.dry, space.drz);
    Log("dk: %f %f %f\n", space.dkx, space.dky, space.dkz);

    int N = space.Nx * space.Ny * space.Nz;

    Log("%ld grid cells.\n", N);


    fftw_init_threads();
    fftw_plan_with_nthreads(nthreads);

    state.psi = (fftw_complex*) fftw_malloc(sizeof(*state.psi) * N);
    state.phi = (fftw_complex*) fftw_malloc(sizeof(*state.phi) * N);
    state.rho = (fftw_complex*) fftw_malloc(sizeof(*state.rho) * N);

    Log("Creating plans...\n");
#define PLAN(sp, st, var, dir) \
    fftw_plan_dft_3d(sp.Nx, sp.Ny, sp.Nz, st. var, st. var, dir, FFTW_ESTIMATE)

    plans.psi_f = PLAN(space, state, psi, FFTW_FORWARD);
    plans.phi_f = PLAN(space, state, phi, FFTW_FORWARD);
    plans.rho_f = PLAN(space, state, rho, FFTW_FORWARD);
    plans.psi_b = PLAN(space, state, psi, FFTW_BACKWARD);
    plans.phi_b = PLAN(space, state, phi, FFTW_BACKWARD);
    plans.rho_b = PLAN(space, state, rho, FFTW_BACKWARD);

    Log("Running simulation...\n");

    ic_test1D(&state, &space);

    write_state(&state, &space);

    for (t=0, step=0; t < tmax; t += dt, step++)
    {
        Log("Time %8.4f  Step %i\n", t, step);

        //----------------------------------------------------------------------
        // Drift
        //----------------------------------------------------------------------
        fftw_execute(plans.psi_f);
        for (i=0; i < N; i++)
        {
            double k2, x,y,z;

            kxyz(i, &space, &x,&y,&z);
            k2 = x*x + y*y + z*z;

            k2 /= 2;

            double c = cos(k2 * dt);
            double s = sin(k2 * dt);
            double p0 = state.psi[i][0];
            double p1 = state.psi[i][1];

            state.psi[i][0] = c*p0 - s*p1;
            state.psi[i][1] = c*p1 + s*p0;
        }
        fftw_execute(plans.psi_b);
        for (i=0; i < N; i++)
        {
            state.psi[i][0] /= N;
            state.psi[i][1] /= N;
        }

        //----------------------------------------------------------------------
        // Gravity
        //----------------------------------------------------------------------
        double prob_tot = 0;
        for (i=0; i < N; i++)
        {
            state.rho[i][0] = pow(hypot(state.psi[i][0], state.psi[i][1]),2);
            state.rho[i][1] = 0;

            state.phi[i][0] = 0;
            state.phi[i][1] = 0;
            prob_tot += state.rho[i][0];
        }
#if 1

        Log("Total Probability is %f\n", prob_tot);
        fftw_execute(plans.rho_f);
        for (i=0; i < N; i++)
        {
            double k2, x,y,z;

            kxyz(i, &space, &x, &y, &z);
            k2 = x*x + y*y + z*z;
            //if (fdim(k2,1e-8) > 0)
            if (fabs(k2) > 1e-12)
            {
                state.phi[i][0] = state.rho[i][0] / k2;
                state.phi[i][1] = 0;
            }
        }
        fftw_execute(plans.phi_b);
        for (i=0; i < N; i++)
        {
            state.phi[i][0] /= N;
            state.phi[i][1] /= N;
        }
#endif

        phi_harmonic_oscillator(&state, &space, 0,0,0, space.xmax/2, 0,0);

        //----------------------------------------------------------------------
        // Kick
        //----------------------------------------------------------------------
        for (i=0; i < N; i++)
        {
            double p0 = state.psi[i][0];
            double p1 = state.psi[i][1];
            double c = cos(state.phi[i][0] * dt);
            double s = sin(state.phi[i][0] * dt);

            state.psi[i][0] = c*p0 - s*p1;
            state.psi[i][1] = c*p1 + s*p0;

            //fprintf(stderr, "p0/p1  %f %f -> %f %f / %f\n", p0, p1, state.psi[i][0], state.psi[i][1], state.phi[i][0]*dt);
        }

#if 0
        fftw_execute(plans.rho_f);
        for (i=0; i < N; i++)
        {
            state.rho[i][0] /= N;
            state.rho[i][1] /= N;
        }
#endif

        write_state(&state, &space);
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

    return EXIT_SUCCESS;
}

