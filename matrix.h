#ifndef MATRIX_H
#define MATRIX_H

#include <stdlib.h>

struct matrix
{
    size_t nr, nc;
    double *v;
    double **r;
};


int alloc_matrix(struct matrix *m, size_t nr, size_t nc);
int write_matrix(struct matrix *m, const char *fmt, ...);

#endif

