#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "matrix.h"

int alloc_matrix(struct matrix *m, size_t nr, size_t nc)
{
    size_t i;

    m->nr = nr;
    m->nc = nc;
    m->v = malloc(sizeof(*m->v) * nr * nc); if (m->v == NULL) return 1;
    m->r = malloc(sizeof(*m->r) * nr); if (m->r == NULL) return 1;

    for (i=0; i < nr; i++)
    {
        m->r[i] = m->v + i*nc;
    }

    return 0;
}

static int write_matrix0(FILE *fp, struct matrix *m)
{
    size_t i,j;

    for (i=0; i < m->nr; i++)
    {
        for (j=0; j < m->nc; j++)
        {
            fprintf(fp, "%24.15e", m->r[i][j]);
        }
        fprintf(fp, "\n");
    }

    return 0;
}

int write_matrix(struct matrix *m, const char *fmt, ...)
{
    char *fname;

    va_list ap;

    va_start(ap, fmt);
    vasprintf(&fname, fmt, ap);
    va_end(ap);

    if (fname == NULL) return 1;

    va_start(ap, fmt);
    vsprintf(fname, fmt, ap);
    va_end(ap);

    FILE *fp = fopen(fname, "wb");
    if (fp != NULL)
    {
        Log("Writing %s...\n", fname);
        write_matrix0(fp, m);
        fclose(fp);
    }

    free(fname);

    return 0;
}
