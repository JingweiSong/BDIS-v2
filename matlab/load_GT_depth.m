function [depth_GT] = load_GT_depth(frameindex,filepath)
%LOAD_GT_DEPTH 此处显示有关此函数的摘要
%   此处显示详细说明


addpath(genpath('pinholeModel'));
img = imread([filepath 'depth_' num2str(frameindex,'%07d') '.png']);
r=double(img(:,:,1));g=double(img(:,:,2));b=double(img(:,:,3));
depth_GT = (b+1)/256+100*(g+1)/256+(r+1)/256;
% depthimage = mediannan(depthimage, wsize);


end

