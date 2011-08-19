#define _GNU_SOURCE
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "becon.h"
#include "log.h"
#include "io_arrays.h"

int write_ascii_rho(struct env *env, const char *fmt, ...)
{
    size_t idx;
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

        fprintf(fp, "%ld %ld %ld", env->state.N, 0L, 0L);
        for (idx=0; idx < env->state.N; idx++)
        {
            double q = pow(env->state.psi[idx][0],2)
                     + pow(env->state.psi[idx][1],2);
            fprintf(fp, "%15.8g\n", q);
        }

        fclose(fp);
    }

    free(fname);

    return 0;
}

