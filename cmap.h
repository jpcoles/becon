#ifndef CMAP_H
#define CMAP_H

typedef void (*cmap_t)(double v, int *r, int *g, int *b);

void cmap_tipsy(double v, int *r, int *g, int *b);
void cmap_astro(double v, int *r, int *g, int *b);

#endif

