function [ data ] = loadImages( parameter_Settings,i)

addpath(genpath(fullfile(parameter_Settings.dataset_path,parameter_Settings.sequence_name)));


    filename = ['depth',int2str(i),'.mat'];
    load(filename);

    data.depthSeries = double(depthimage);

    
    
    filename = ['RGB_image',int2str(i),'.mat'];
    if(exist(filename,'file')==0)
        data.RGB_image = zeros([size(data.depthSeries),3],'uint8');
    else
        load(filename);
        data.RGB_image = RGB_image;
    end

end

