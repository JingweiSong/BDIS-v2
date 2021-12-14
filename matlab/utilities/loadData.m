function [ data ] = loadData( sequence_name, depth_datapath,start_frame, end_frame )


if ~exist('sequence_name','var')
    error('No such directory');
end


addpath(genpath(fullfile(depth_datapath,sequence_name)));

% fileName = 'left.mp4'; 
% obj = VideoReader(fileName);

for i = start_frame : end_frame
    filename = ['depth',int2str(i),'.mat'];
    load(filename);
    data.depthSeries{i} = depthimage;

%     filename = ['depthinfo',int2str(i),'.mat'];
%     load(filename);
%     data.depthinfo{i} = datainfo_depth;
    filename = ['weightmap',int2str(i),'.mat'];
    load(filename);
    data.weightmap{i} = weightmap;
    
    filename = ['RGB_image',int2str(i),'.mat'];
    if(exist(filename,'file')==0)
        data.RGB_image{i} = zeros([size(data.depthSeries{i}),3],'uint8');
    else
        load(filename);
        data.RGB_image{i} = RGB_image;
        
%         data.RGB_image{i} = read(obj,i);
    end
end

end

