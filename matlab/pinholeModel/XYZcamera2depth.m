function [ u,v,value ] = XYZcamera2depth( K,pointcloud )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script:   Transform a pointcloudcloud [X,Y,Z] to depth image [u,v,value]
%   Method:
%   Input:
%             K:                Camera intrinsic matrix
%             pointcloud:            1¡Á3 pointcloud
%   Returns:
%             u:                the 'u'th row of image
%             v:                the 'v'th column of image
%             value:            value of depth image [u,v]
%   Author:   Jingwei Song.     02/08/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_point = size(pointcloud,1);
v       = zeros(num_point,1);
u       = zeros(num_point,1);
value   = zeros(num_point,1);

for i = 1 : num_point
    v(i)      = round(K(1,1)*(pointcloud(i,1)./pointcloud(i,3)) + K(1,3));
    u(i)      = round(K(2,2)*(pointcloud(i,2)./pointcloud(i,3)) + K(2,3));
    value(i)  = pointcloud(i,3);
end
end

