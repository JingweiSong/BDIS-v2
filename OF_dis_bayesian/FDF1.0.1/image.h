#ifndef __IMAGE_H_
#define __IMAGE_H_

#include <stdio.h>

#define MIN_TA(a, b) ((a) < (b) ? (a) : (b))
#define MAX_TA(a, b) ((a) > (b) ? (a) : (b))
#define MINMAX_TA(a,b) MIN_TA( MAX_TA(a,0) , b-1 )

#ifdef __cplusplus
extern "C" {
#endif





typedef struct image_s
{
  int width;
  int height;
  int stride;
  float *c1;
} image_t;


typedef struct color_image_s
{
    int width;
    int height;
    int stride;
    float *c1;
    float *c2;
    float *c3;
} color_image_t;


typedef struct color_image_pyramid_s 
{
  float scale_factor;
  int min_size;
  int size;
  color_image_t **images;
} color_image_pyramid_t;


typedef struct convolution_s
{
    int order;
    float *coeffs;
    float *coeffs_accu;
} convolution_t;




image_t *image_new(const int width, const int height);


image_t *image_cpy(const image_t *src);


void image_erase(image_t *image);


void image_delete(image_t *image);


void image_mul_scalar(image_t *image, const float scalar);


color_image_t *color_image_new(const int width, const int height);


color_image_t *color_image_cpy(const color_image_t *src);


void color_image_erase(color_image_t *image);


void color_image_delete(color_image_t *image);


void resize_if_needed_newsize(image_t *im, const int w, const int h);




image_t *image_resize_bilinear(const image_t *src, const float scale);


void image_resize_bilinear_newsize(image_t *dst, const image_t *src, const int new_width, const int new_height);


color_image_t *color_image_resize_bilinear(const color_image_t *src, const float scale);




float *gaussian_filter(const float sigma, int *fSize);


convolution_t *convolution_new(int order, const float *half_coeffs, const int even);


void convolve_horiz(image_t *dest, const image_t *src, const convolution_t *conv);


void convolve_vert(image_t *dest, const image_t *src, const convolution_t *conv);


void convolution_delete(convolution_t *conv);


void color_image_convolve_hv(color_image_t *dst, const color_image_t *src, const convolution_t *horiz_conv, const convolution_t *vert_conv);


void image_convolve_hv(image_t *dst, const image_t *src, const convolution_t *horiz_conv, const convolution_t *vert_conv);





color_image_pyramid_t *color_image_pyramid_create(const color_image_t *src, const float scale_factor, const int min_size, const float spyr);


void color_image_pyramid_delete(color_image_pyramid_t *pyr);

#ifdef __cplusplus
}
#endif


#endif


