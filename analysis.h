#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "becon.h"
#include "matrix.h"
#include "analysis.h"

int radial_avg_density(struct env *env, struct matrix *m, double x, double y, double z);
int radial_avg_potential(struct env *env, struct matrix *m, double x, double y, double z);

#endif
