#ifndef IO_IMAGE_H
#define IO_IMAGE_H

#include "frame_buffer.h"

int write_image(const char *fmt, int step, struct frame_buffer *fb);
int capture_image(struct space *space, struct state *state, struct frame_buffer *fb);

#endif
