function [ pointcloud_new ] = getObservablePointcloud( pointcloud,voxelsize,K )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script:   Get observable point cloud
%   Method:   using pin hole camera to predict observable pointcloud
%   Input:    
%             pointcloud:       3¡ÁN pointcloud
%             voxelsize:        2¡Á1 Voxel size(like kinect 480¡Á240)
%             voxelunit:        Size of a voxel. 
%             K:                Camera intrinsic matrix
%   Returns:  
%             pointcloud_new:   Observable point cloud
%             faces_new:        New Triangle faces
%   Author:   Jingwei Song.   02/08/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% observable_range(1,1) = - voxelsize(1) * voxelunit / 2;
% observable_range(1,2) = observable_range(1,1) + (voxelsize(1)-1) * voxelunit;
% observable_range(2,1) = - voxelsize(2) * voxelunit / 2;
% observable_range(2,2) = observable_range(2,1) + (voxelsize(2)-1) * voxelunit;

px = round(K(1,1)*(pointcloud(:,1)./pointcloud(:,3)) + K(1,3));
py = round(K(2,2)*(pointcloud(:,2)./pointcloud(:,3)) + K(2,3));
isValid = (1<=px & px <= voxelsize(1) & 1<=py & py<= voxelsize(2));
pointcloud_new = pointcloud(isValid,:);
end

