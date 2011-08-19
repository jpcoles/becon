#ifndef IO_TIPSY_H
#define IO_TIPSY_H

#include <stdint.h>
#include "becon.h"

struct tipsy_header
{
    double time;
    uint32_t nBodies;
    uint32_t nDims;
    uint32_t nGas;
    uint32_t nDark;
    uint32_t nStar;
};

struct tipsy_dm
{
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
};


int write_tipsy_grid(struct env *env, const char *fmt, ...);

#endif

