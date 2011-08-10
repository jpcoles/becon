#ifndef IO_GRAFIC_H
#define IO_GRAFIC_H

#include "becon.h"

int alloc_grafic(const char *dirname, struct state *state, struct space *space);
int read_grafic(const char *dirname, struct state *state, struct space *space);

#endif
