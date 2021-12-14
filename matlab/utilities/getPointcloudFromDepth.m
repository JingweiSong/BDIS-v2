function [ currentdepth_pts ] = getPointcloudFromDepth( cameraIntrinsicParam,depth,RGB_image )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: Get point cloud from depth
%   Method:   As the depth image comes from pin-hole camera, we first
%            transform it to CMAERA coordinate then to world cordinate.
%   Input:    cameraIntrinsicParam:      Camera intrinsic matrix
%             depth:          Orthographical depth image
%                             voxel
%             RGB_image:      3 channels RGB image
%             volume_limit:   record limits for points in a volume; 
%                             COL1 minXYZ; COL2 maxXYZ
%   Returns:  currentdepth_pts: Pointcloud of current depth
%   Author:   Jingwei Song.   23/08/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[num_imagerow,num_imagecolunm] = size(depth);
[currentdepth_pts] = RGBD2XYZcamera(cameraIntrinsicParam,depth,RGB_image,num_imagerow,num_imagecolunm);


end

