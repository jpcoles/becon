#ifndef IO_IMAGE_H
#define IO_IMAGE_H

#include "frame_buffer.h"
#include "cmap.h"

int write_image(struct frame_buffer *fb, const char *fmt, ...);
int capture_image(struct env *env, struct frame_buffer *fb);
int capture_image_log(double min, double max, cmap_t cmap, struct env *env, struct frame_buffer *fb);

#endif
