#include <assert.h>
#include <stdlib.h>
#include "frame_buffer.h"

int alloc_frame_buffer(struct frame_buffer *fb, int32_t width, int32_t height)
{
    fb->samples = 3;
    fb->width  = width;
    fb->height = height;
    fb->size = width * height;

    fb->buf = malloc(fb->size * sizeof(double));
    fb->img = malloc(fb->size * fb->samples);

    if (fb->img == NULL) return 1;
    if (fb->buf == NULL) return 1;

    return 0;
}

int free_frame_buffer(struct frame_buffer *fb)
{
    assert(fb->buf != NULL);
    assert(fb->img != NULL);
    free(fb->buf);
    free(fb->img);

    return 0;
}

void write_frame_buffer(FILE *outfile, struct frame_buffer *fb)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width      = fb->width;
    cinfo.image_height     = fb->height;
    cinfo.input_components = fb->samples;
    cinfo.in_color_space   = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_start_compress(&cinfo, TRUE);

    JSAMPROW row_pointer[1];              /* pointer to a single row */
    int row_stride;                       /* physical row width in buffer */

    row_stride = cinfo.image_width * fb->samples;   /* JSAMPLEs per row in image_buffer */

    while (cinfo.next_scanline < cinfo.image_height) {
        row_pointer[0] = & (fb->img[(cinfo.image_height-cinfo.next_scanline-1) * row_stride]);
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

}
