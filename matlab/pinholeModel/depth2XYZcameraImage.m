function [XYZcamera,XYZnormal,d_depth_X,d_depth_Y] = depth2XYZcameraImage(K, depth,image_size)
    XYZcamera = zeros(image_size(1),image_size(2),3);
    XYZnormal = zeros(image_size(1),image_size(2),3);

    [x,y] = meshgrid(1:image_size(2), 1:image_size(1));
    XYZcamera(:,:,1) = x;
    XYZcamera(:,:,2) = y;
    XYZcamera(:,:,3) = nan;
    
    %   Estimating normals
    [ d_depth_X,d_depth_Y] = calculateDerivativeImage(depth,'robert');
    XYZnormal(:,:,1) = -1 * d_depth_X;
    XYZnormal(:,:,2) = -1 * d_depth_Y;
    XYZnormal(:,:,3) = -1;
    %   Regularize normal
    magnitude = sqrt(XYZnormal(:,:,1).^2 + XYZnormal(:,:,2).^2 + XYZnormal(:,:,3).^2);
    XYZnormal(:,:,1) = XYZnormal(:,:,1) ./ magnitude;
    XYZnormal(:,:,2) = XYZnormal(:,:,2) ./ magnitude;
    XYZnormal(:,:,3) = XYZnormal(:,:,3) ./ magnitude;
    
    %firstly find the depth==0;
    [rows,cols] = find(depth~=0 & ~isnan(depth) & ~isinf(depth)); 
    index_depth_image = sub2ind(image_size,rows,cols);
    index_xyzcamera_x = sub2ind([image_size 3],rows,cols,ones(size(rows,1),1)*1);
    index_xyzcamera_y = sub2ind([image_size 3],rows,cols,ones(size(rows,1),1)*2);
    index_xyzcamera_z = sub2ind([image_size 3],rows,cols,ones(size(rows,1),1)*3);
    XYZcamera(index_xyzcamera_x) = (x(index_depth_image)-K(1,3)).*depth(index_depth_image)/K(1,1);
    XYZcamera(index_xyzcamera_y) = (y(index_depth_image)-K(2,3)).*depth(index_depth_image)/K(2,2);
    XYZcamera(index_xyzcamera_z) = depth(index_depth_image);
    XYZcamera(:,:,4) = depth~=0;
end
