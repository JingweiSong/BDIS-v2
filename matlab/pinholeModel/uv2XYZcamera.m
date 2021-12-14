function [ XYZcamera ] = uv2XYZcamera( K, u,v,value_uv,image_row,image_col )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: Convert (u,v) of a depth image to a single point
%   Method:   Similar to depth2XYZcamera
%   Input:    K:                Intrinsic parameter
%             u,v:              Pixel location
%             image_row:        Row of image
%             image_col:        Colunmn of image
%   Returns:  XYZcamera:        Converted image
%   Author:   Jingwei Song.     19/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Throw exception
if((u>image_row)||(u<0)||(u>image_row)||(u<0))
    error('(u,v) not legal');
end
XYZcamera = zeros(1,3);
XYZcamera(1) = (v-K(1,3)).*value_uv/K(1,1);
XYZcamera(2) = (u-K(2,3)).*value_uv/K(2,2);
XYZcamera(3) = value_uv;

end

