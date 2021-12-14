function [ pointcloud_w,pose_matrix ] = CameraCord2world( pointcloud_c,camera_center_w,deviation_angle )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script:   (i) Convert camera_center center to [0,0,-sqrt(x^2+y^2+z^2)]
%                   x,y,z defined by camera_center_w
%             (ii) Calcuate rotation R between(current camera_center to
%             camera_center_w]
%             (iii)Rotate the model with R 
%   Input:    
%             pointcloud_c:     3¡ÁN pointcloud (camera cordinate)
%             camera_center_w:  camera center (world cordinate)
%             deviation_angle:  camera to masscenter deviation angle.
%   Returns:  
%             pointcloud_w:     3¡ÁN pointcloud (world cordinate)
%             pose_matrix:      Transformation matrix in the form of [R T;0 1];
%   Author:   Jingwei Song.   17/08/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
camera_center_tmp = [0,0, ...
                    -sqrt(camera_center_w(1)^2 + camera_center_w(2)^2 ...
                    + camera_center_w(3)^2)];
pointcloud_w = pointcloud_c + repmat(camera_center_tmp,[size(pointcloud_c,1),1]);
%camera_center_tmp =[0,0,0] - center_mass_c;


r               = vrrotvec(camera_center_tmp,camera_center_w);
rotate_matrix   = vrrotvec2mat(r);

%   Default observation direction is from camera center to model center
%   alpha beta gamma defines minor deviations
deviation_rot   = eul2rotm(deviation_angle);

pointcloud_w    = deviation_rot*rotate_matrix * pointcloud_w';
pointcloud_w    = pointcloud_w';

pose_matrix     = zeros(4,4);
pose_matrix(1:3,1:3) = deviation_rot * rotate_matrix;
pose_matrix(1:3,end) = (pose_matrix(1:3,1:3)*camera_center_tmp')';
pose_matrix(4,  4)   = 1;

end

