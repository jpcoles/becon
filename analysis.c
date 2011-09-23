#include "becon.h"
#include "log.h"
#include "matrix.h"
#include "analysis.h"

int radial_avg_density(struct env *env, struct matrix *m, double x, double y, double z)
{
    int32_t i,j,k;
    size_t idx=0;
    struct space *sp = &env->space;

    if (m->nr != sp->Nmax/2+1 || m->nc != 3)
    {
        Log("WARNING: radial_avg_density() : Matrix dimensions %ldx%ld should be %ldx%ld.\n", m->nr, m->nc, sp->Nmax/2+1, 3);
        Log("WARNING: radial_avg_density() : Density not calculated.\n");
        return 1;
    }

    for (idx=0; idx < m->nr; idx++)
    {
        m->r[idx][0] = idx;
        m->r[idx][1] = 0;
        m->r[idx][2] = 0;
    }

    idx = 0;
    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->dx + sp->xmin;
        double y0 = j*sp->dx + sp->ymin;
        double z0 = k*sp->dx + sp->zmin;

        int r = sqrt(pow(x0,2) + pow(y0,2) + pow(z0,2)) / sp->dx;

        if (r <= sp->Nmax/2)
        {
            m->r[r][1] += sqrt(cmag2(env->state.psi[idx]));
            m->r[r][2]++;
        }

        idx++;
    }

    for (idx=0; idx < m->nr; idx++)
    {
        if (m->r[idx][2] != 0)
            m->r[idx][1] /= m->r[idx][2];
    }

    return 0;
}

int radial_avg_potential(struct env *env, struct matrix *m, double x, double y, double z)
{
    int32_t i,j,k;
    size_t idx=0;
    struct space *sp = &env->space;

    if (m->nr != sp->Nmax/2+1 || m->nc != 3)
    {
        Log("WARNING: radial_avg_potential() : Matrix dimensions %ldx%ld should be %ldx%ld.\n", m->nr, m->nc, sp->Nmax/2+1, 3);
        Log("WARNING: radial_avg_potential() : Potential not calculated.\n");
        return 1;
    }

    for (idx=0; idx < m->nr; idx++)
    {
        m->r[idx][0] = idx;
        m->r[idx][1] = 0;
        m->r[idx][2] = 0;
    }

    idx = 0;
    for (k=0; k < sp->Nz; k++)
    for (j=0; j < sp->Ny; j++)
    for (i=0; i < sp->Nx; i++)
    {
        double x0 = i*sp->dx + sp->xmin;
        double y0 = j*sp->dx + sp->ymin;
        double z0 = k*sp->dx + sp->zmin;

        int r = sqrt(pow(x0,2) + pow(y0,2) + pow(z0,2)) / sp->dx;

        if (r <= sp->Nmax/2)
        {
            m->r[r][1] += env->state.phi[idx][0];
            m->r[r][2]++;
        }

        idx++;
    }

    for (idx=0; idx < m->nr; idx++)
    {
        if (m->r[idx][2] != 0)
            m->r[idx][1] /= m->r[idx][2];
    }

    return 0;
}
