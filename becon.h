#ifndef BECON_H
#define BECON_H

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

struct options
{
    int32_t with_gravity;
    char *input;
    double dt;
    double tmax;
};

#endif
