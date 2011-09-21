#ifndef BECON_H
#define BECON_H

#include <sys/types.h>
#include <fftw3.h>
#include <math.h>

struct state
{
    size_t N, vN;
    fftw_complex *psi;
    fftw_complex *phi;
    fftw_complex *rho;
    fftw_complex *S;
    fftw_complex *vphi;
};

struct space
{
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    int32_t Nx, Ny, Nz;
    int32_t vNx, vNy, vNz;    /* For vacuum conditions */
    int32_t Nmax;
    double dx;
    double dk;
};

struct plans
{
    fftw_plan psi_f;
    fftw_plan phi_f;
    fftw_plan rho_f;
    fftw_plan S_f;

    fftw_plan psi_b;
    fftw_plan phi_b;
    fftw_plan rho_b;
    fftw_plan vphi_b;
};

struct options
{
    int32_t with_gravity;
    int32_t with_vacuum;
    char *input;
    double dt;
    double tmax;
    int noutputs;
};

struct cosmo
{
    double a;
    double a_start;
    double H0;
    double h;

    double rho_crit;

    /* density parameters */
    double omega_m, omega_v, omega_r, omega_k;

    /* hbar / m */
    //double h_m;
    //double h,m;


};

struct bec
{
    double m;
};

struct consts
{
    struct {
        double c;
        double hbar;
        double H0;
        double G;
    } si;

    struct {
        double c;
        double hbar;
        double H0;
        double G;
    } in;
};

struct env
{
    struct space space;
    struct state state;

    struct cosmo cosmo;

    struct bec bec;
    struct consts consts;

    struct options opts;

    double rho_max, rho_min, rho_avg;
    double phi_max;

    double kick_exp, drift_exp;
    double eta;
    double dt;

    double kmax2;
};

static inline int imin(int a, int b)
{
    if (a < b) return a;
    return b;
}

static inline int imax(int a, int b)
{
    if (a > b) return a;
    return b;
}

static inline double MIN(double a, double b)
{
    if (a < b) return a;
    return b;
}

static inline double MAX(double a, double b)
{
    if (a > b) return a;
    return b;
}

static inline double cmag2(fftw_complex a)
{
    return pow(a[0],2) + pow(a[1],2);
}

#endif
