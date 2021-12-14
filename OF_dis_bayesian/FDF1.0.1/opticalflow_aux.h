#ifndef __OPTICALFLOW_AUX_
#define __OPTICALFLOW_AUX_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

#include "image.h"

#if (SELECTCHANNEL==1 | SELECTCHANNEL==2)

void image_warp(image_t *dst, image_t *mask, const image_t *src, const image_t *wx, const image_t *wy);


void get_derivatives(const image_t *im1, const image_t *im2, const convolution_t *deriv, image_t *dx, image_t *dy, image_t *dt, image_t *dxx, image_t *dxy, image_t *dyy, image_t *dxt, image_t *dyt);

#else

void image_warp(color_image_t *dst, image_t *mask, const color_image_t *src, const image_t *wx, const image_t *wy);


void get_derivatives(const color_image_t *im1, const color_image_t *im2, const convolution_t *deriv, color_image_t *dx, color_image_t *dy, color_image_t *dt, color_image_t *dxx, color_image_t *dxy, color_image_t *dyy, color_image_t *dxt, color_image_t *dyt);

#endif
    
    

void compute_smoothness(image_t *dst_horiz, image_t *dst_vert, const image_t *uu, const image_t *vv, const convolution_t *deriv_flow, const float quarter_alpha);



void sub_laplacian(image_t *dst, const image_t *src, const image_t *weight_horiz, const image_t *weight_vert);


void compute_data_and_match(image_t *a11, image_t *a12, image_t *a22, image_t *b1, image_t *b2, image_t *mask, image_t *wx, image_t *wy, image_t *du, image_t *dv, image_t *uu, image_t *vv, color_image_t *Ix, color_image_t *Iy, color_image_t *Iz, color_image_t *Ixx, color_image_t *Ixy, color_image_t *Iyy, color_image_t *Ixz, color_image_t *Iyz, image_t *desc_weight, image_t *desc_flow_x, image_t *desc_flow_y, const float half_delta_over3, const float half_beta, const float half_gamma_over3);


#if (SELECTCHANNEL==1 | SELECTCHANNEL==2)
void compute_data(image_t *a11, image_t *a12, image_t *a22, image_t *b1, image_t *b2, image_t *mask, image_t *wx, image_t *wy, image_t *du, image_t *dv, image_t *uu, image_t *vv, image_t *Ix, image_t *Iy, image_t *Iz, image_t *Ixx, image_t *Ixy, image_t *Iyy, image_t *Ixz, image_t *Iyz, const float half_delta_over3, const float half_beta, const float half_gamma_over3);
#else
void compute_data(image_t *a11, image_t *a12, image_t *a22, image_t *b1, image_t *b2, image_t *mask, image_t *wx, image_t *wy, image_t *du, image_t *dv, image_t *uu, image_t *vv, color_image_t *Ix, color_image_t *Iy, color_image_t *Iz, color_image_t *Ixx, color_image_t *Ixy, color_image_t *Iyy, color_image_t *Ixz, color_image_t *Iyz, const float half_delta_over3, const float half_beta, const float half_gamma_over3);
#endif

#if (SELECTCHANNEL==1 | SELECTCHANNEL==2)
void compute_data_DE(image_t *a11, image_t *b1, image_t *mask, image_t *wx, image_t *du, image_t *uu, image_t *Ix, image_t *Iy, image_t *Iz, image_t *Ixx, image_t *Ixy, image_t *Iyy, image_t *Ixz, image_t *Iyz, const float half_delta_over3, const float half_beta, const float half_gamma_over3);
#else
void compute_data_DE(image_t *a11, image_t *b1, image_t *mask, image_t *wx, image_t *du, image_t *uu, color_image_t *Ix, color_image_t *Iy, color_image_t *Iz, color_image_t *Ixx, color_image_t *Ixy, color_image_t *Iyy, color_image_t *Ixz, color_image_t *Iyz, const float half_delta_over3, const float half_beta, const float half_gamma_over3);
#endif



void descflow_resize(image_t *dst_flow_x, image_t *dst_flow_y, image_t *dst_weight, const image_t *src_flow_x, const image_t *src_flow_y, const image_t *src_weight);


void descflow_resize_nn(image_t *dst_flow_x, image_t *dst_flow_y, image_t *dst_weight, const image_t *src_flow_x, const image_t *src_flow_y, const image_t *src_weight);

#ifdef __cplusplus
}
#endif


#endif
