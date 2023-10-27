close all;clear all;clc;

addpath(genpath('utilities'));
addpath(genpath('pinholeModel'));
parameter_Settings.dataset_path     ='../dataset/';


Dataindex = 30;
switch Dataindex
    case 30
        parameter_Settings.end_frame = 10;
    case 31
        parameter_Settings.end_frame = 15;
    case 32
        parameter_Settings.end_frame = 25;
    case 33
        parameter_Settings.end_frame = 25;
    otherwise
        disp('Wrong index')
end
parameter_Settings.GTdepth_path     = [parameter_Settings.dataset_path num2str(Dataindex) '\'];
parameter_Settings.sequence_name    = [num2str(Dataindex) '\points'];
parameter_Settings.camera_intrinsic  = [520	0	320;...
    0	520	240;...
    0	0	1];
parameter_Settings.base_line = 0.5;
parameter_Settings.start_frame =1;
frame_gap = 1;

 
frameindex = parameter_Settings.start_frame;                        
while ( frameindex <= parameter_Settings.end_frame  )
%     close all;

    result_Bayesian = double(imread([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/bayesian_output' int2str(frameindex) '.png']))/10;
    

    Data.RGB_image = imread([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/left/' int2str(frameindex) '.png']);
    
    %   Filter invalid pixels
    result_Bayesian(Data.RGB_image(:,:,1)==0) = NaN;


    image_size = size(Data.RGB_image);image_size=image_size(:,1:2);

    depth_Bayesian = getDepthFromDisparitymap(result_Bayesian,parameter_Settings,image_size);
    depth_Bayesian(depth_Bayesian>30) = NaN;
    [currentdepth_pts_Bayesian] = getPointcloudFromDepth( ...
                                parameter_Settings.camera_intrinsic,...
                                depth_Bayesian,...
                                Data.RGB_image);
                            
    figure('Name',['Left image' int2str(frameindex)],'NumberTitle','off');
    imshow(Data.RGB_image);
    figure('Name',['BDIS shape' int2str(frameindex)],'NumberTitle','off');
    pcshow(currentdepth_pts_Bayesian(:,1:3),uint8(currentdepth_pts_Bayesian(:,4:6)),'MarkerSize',20);set(gca,'visible','off'); set(gcf,'color','w');
  
    
    if(exist([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/prob_out_output' int2str(frameindex)],'file'))
        result_prob = importdata([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/prob_out_output' int2str(frameindex)]);
        depthGT = load_GT_depth(frameindex,[parameter_Settings.dataset_path num2str(Dataindex) '/']);
        error = abs(depthGT-depth_Bayesian);
        bool_ind = (result_prob~=0 ) & (error<10);
        num_valid_all = sum(sum(bool_ind));
        error = error(bool_ind);
        result_prob = result_prob(bool_ind);
        f = parameter_Settings.camera_intrinsic(1,1);
        base_line = parameter_Settings.base_line;
        var_result_Bayesian = (f*base_line)^2*result_prob ./ result_Bayesian(bool_ind).^4;
        num_valid = sum(error<1.96*sqrt(var_result_Bayesian));
        rate_1_96sigma = num_valid/num_valid_all;
        disp(['The average coverage rate (1.96 sigma bound):  ' num2str(rate_1_96sigma)]);
    end
    frameindex = frameindex + frame_gap
end

function [depth_GT] = load_GT_depth(frameindex,filepath)

addpath(genpath('pinholeModel'));
img = imread([filepath 'depth_' num2str(frameindex,'%07d') '.png']);
r=double(img(:,:,1));g=double(img(:,:,2));b=double(img(:,:,3));
depth_GT = (b+1)/256+100*(g+1)/256+(r+1)/256- 0.5;
end
