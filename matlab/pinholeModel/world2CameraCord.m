function [ pointcloud_c,pose_matrix ] = world2CameraCord( pointcloud_w,camera_center_w,deviation_angle )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script:   (i) Convert camera center or original point
%             (ii) Calcuate rotation R between(center of mass[x,y,z] to
%             [0,0,sqrt(x^2+y^2+z^2)]
%             (iii)Rotate the model with R 
%             (iv) Rotate the model with deviation_matrix
%   Input:    
%             pointcloud_w:     3¡ÁN pointcloud (world cordinate)
%             camera_center_w:  camera center (world cordinate)
%             deviation_angle:  camera to masscenter deviation angle.
%   Returns:  
%             pointcloud_c:     3¡ÁN pointcloud (world cordinate)
%   Author:   Jingwei Song.   17/08/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pointcloud_c = pointcloud_w - repmat(camera_center_w,[size(pointcloud_w,1),1]);
center_mass_source  = - camera_center_w;
center_mass_target  = [0,0,sqrt(camera_center_w(1)^2 + ...
                                camera_center_w(2)^2 + camera_center_w(3)^2)];


r               = vrrotvec(center_mass_source,center_mass_target);
rotate_matrix   = vrrotvec2mat(r);

pointcloud_c    = rotate_matrix * pointcloud_c';

%   Default observation direction is from camera center to model center
%   alpha beta gamma defines minor deviations
deviation_rot   = eul2rotm(deviation_angle);

pointcloud_c     = deviation_rot * pointcloud_c;

pointcloud_c    = pointcloud_c';

pose_matrix     = zeros(4,4);
pose_matrix(1:3,1:3) = rotate_matrix * deviation_rot;
pose_matrix(1:3,end) = (pose_matrix(1:3,1:3)*(-camera_center_w)')';
pose_matrix(4,  4)   = 1;
end

