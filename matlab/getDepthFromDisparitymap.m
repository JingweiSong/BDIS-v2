function [ depth_image ] = getDepthFromDisparitymap( OF_image,parameter_Settings,image_size )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% OF_image = OF_image';
% OF_image = resize_disparity(OF_image,image_size);
f = parameter_Settings.camera_intrinsic(1,1);
depth_image = zeros(image_size);
ind = OF_image~=0;
depth_image(ind) = parameter_Settings.base_line*f./OF_image(ind);
depth_image(~ind) = nan;
end

function disparity_map = resize_disparity(disparity_map_tmp,size_image)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: scale the disparity map to the image size;
%   Method:   the resize function in matlab work not well.
%   Input:    disparity_map_tmp: the original out put of the elas
%   algorithm;
%   Returns:
%   Author:   Jun Wang.   08/09/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disparity_map = zeros(size_image);
size_disparity = size(disparity_map_tmp);
scales_size = size_disparity./size_image;
[rows_image,cols_image]=meshgrid(1:size_image(1),1:size_image(2));
rows_image = reshape(rows_image,[numel(rows_image),1]);
cols_image = reshape(cols_image,[numel(cols_image),1]);
rows_disparity = round(rows_image.*scales_size(1));
cols_disparity = round(cols_image.*scales_size(2));
is_invalid = rows_disparity<1 | rows_disparity>size_disparity(1) | cols_disparity < 1 | cols_disparity > size_disparity(2);
rows_disparity(is_invalid) = [];
cols_disparity(is_invalid) =[];
rows_image(is_invalid) =[];
cols_image(is_invalid)=[];
index_image = sub2ind(size_image,rows_image,cols_image);
index_disparity = sub2ind(size_disparity,rows_disparity,cols_disparity);

disparity_map(index_image) = disparity_map_tmp(index_disparity);
end