close all;clear all;clc;

addpath(genpath('utilities'));
addpath(genpath('pinholeModel'));
parameter_Settings.dataset_path     ='../dataset/';


Dataindex = 6;
if(Dataindex == 6)
    parameter_Settings.camera_intrinsic  = [6.958775809222639737e+02 0.000000000000000000e+00 2.804026432037353516e+02;0.000000000000000000e+00 6.958775809222639737e+02 2.523805904388427734e+02;0 0  1.000000];
    parameter_Settings.base_line = 5.379;
    parameter_Settings.start_frame = 1;
    parameter_Settings.end_frame = 1051;
    frame_gap = 50;
elseif(Dataindex == 20)
    parameter_Settings.camera_intrinsic = [3.437898110138030461e+02 0.000000000000000000e+00 3.629623012542724609e+02;...
        0.000000000000000000e+00 3.437898110138030461e+02 1.555652084350585938e+02;0 0  1.000000];
    parameter_Settings.base_line = 5.200;
    parameter_Settings.start_frame = 801;
    parameter_Settings.end_frame = 2401;
    frame_gap = 100;
elseif(Dataindex == 21)
    parameter_Settings.camera_intrinsic = [3.437898110138030461e+02 0.000000000000000000e+00 3.629623012542724609e+02;...
        0.000000000000000000e+00 3.437898110138030461e+02 1.555652084350585938e+02;0 0  1.000000];
    parameter_Settings.base_line = 5.200;
    parameter_Settings.start_frame = 1001;
    parameter_Settings.end_frame   = 4201;%4201;
    frame_gap = 100;
end
       

frameindex = parameter_Settings.start_frame;                        
while ( frameindex <= parameter_Settings.end_frame  )
%     close all;
    
   
    
    result_Bayesian = importdata([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/bayesian_output' int2str(frameindex)]);
    
    Data.RGB_image = imread([parameter_Settings.dataset_path num2str(Dataindex) '/images_rectified/left/' int2str(frameindex) '.png']);

    
    %   Filter invalid pixels
    result_Bayesian(Data.RGB_image(:,:,1)==0) = NaN;
   
    
    image_size = size(Data.RGB_image);image_size=image_size(:,1:2);
    depth_Bayesian = getDepthFromDisparitymap(result_Bayesian,parameter_Settings,image_size);
    depth_Bayesian(depth_Bayesian>300) = NaN;
    [currentdepth_pts_Bayesian] = getPointcloudFromDepth( ...
                                parameter_Settings.camera_intrinsic,...
                                depth_Bayesian,...
                                Data.RGB_image);

    figure('Name',['Left image' int2str(frameindex)],'NumberTitle','off');
    imshow(Data.RGB_image);
    figure('Name',['BDIS shape' int2str(frameindex)],'NumberTitle','off');
    pcshow(currentdepth_pts_Bayesian(:,1:3),uint8(currentdepth_pts_Bayesian(:,4:6)),'MarkerSize',20);set(gca,'visible','off'); set(gcf,'color','w');
 
   
    frameindex = frameindex + frame_gap
end
