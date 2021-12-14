function [ depthimage ] = points2depth( pointcloud,face,num_imagerow,num_imagecol,cameraIntrinsicParam )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: Points to depth(orthographical projection) 
%   Method:   See readme
%   Input:    pointcloud:       Current point cloud
%             face:             Faces of model
%             num_imagerow:     Row of image
%             num_imagecol:     Column of image
%   Returns:   
%             depthimage:       Depth image
%             datainfo_depth:   Depth image info including
%                               [num_imagerow, num_imagecol,0;
%                                incremental_x,incremental_y,0;
%                                max_x,       max_y,     max_z;
%                                min_x,       min_y,     min_z;
%                                shift_origin(3),0,   0;] 
%   Author:   Jingwei Song.     25/07/2016
%   Modify:   Jingwei Song.     16/08/2016. Delete rotation. This has been
%                               done in world2CamearCord 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pointcloud = pointcloud';

%===============Test Code===============%
pointcloud_new  = getObservablePointcloud(pointcloud,[num_imagerow,num_imagecol],cameraIntrinsicParam);   
pc_cloud = pointCloud(pointcloud_new);
figure;showPointCloud(pc_cloud);  
clear pc_cloud pointcloud_new;
%=======================================%

depthimage = ProjectModel(face, pointcloud, num_imagerow,num_imagecol,cameraIntrinsicParam);



for i = 1 : num_imagerow
    for j = 1 : num_imagecol
        if(depthimage(i,j) == 0)
            depthimage(i,j) = nan;
        end
    end
end




end

