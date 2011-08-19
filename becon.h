#ifndef BECON_H
#define BECON_H

#include <sys/types.h>
#include <fftw3.h>

struct state
{
    size_t N;
    fftw_complex *psi;
    fftw_complex *phi;
    fftw_complex *rho;
};

struct space
{
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    int32_t Nx, Ny, Nz;
    double dx;
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

struct options
{
    int32_t with_gravity;
    char *input;
    double dt;
    double tmax;
};

struct cosmo
{
    double a_start;
    double H0;

    double rho_crit;

    /* density parameters */
    double omega_m, omega_v, omega_r, omega_k;

    /* hbar / m */
    double h_m;
    double h,m;
};

struct env
{
    struct space space;
    struct state state;

    struct cosmo cosmo;

    struct options opts;

    double rho_max, rho_min;
};

#endif
