function [xyzPoints, xyzNormal ,xyzWeight] = RGBD2XYZcamera(K, depth,RGB_image,image_row,image_col,weight)
% column:1-3:positions 4-6:RGB 7-9:normals 10:weight

    [x_size,y_size]  = size(depth);
    [x,y] = meshgrid(1:y_size, 1:x_size);
    XYZcamera(:,:,1) = (x-K(1,3)).*depth/K(1,1);
    XYZcamera(:,:,2) = (y-K(2,3)).*depth/K(2,2);
    XYZcamera(:,:,3) = depth;
    XYZcamera(:,:,4) = depth~=0;
    
    x = XYZcamera(:,:,1);
    y = XYZcamera(:,:,2);
    z = XYZcamera(:,:,3);
    xyzPoints(:,1)=reshape(x,[image_row*image_col,1]);
    xyzPoints(:,2)=reshape(y,[image_row*image_col,1]);
    xyzPoints(:,3)=reshape(z,[image_row*image_col,1]);
    xyzPoints(:,4)=reshape(RGB_image(:,:,1),[image_row*image_col,1]);
    xyzPoints(:,5)=reshape(RGB_image(:,:,2),[image_row*image_col,1]);
    xyzPoints(:,6)=reshape(RGB_image(:,:,3),[image_row*image_col,1]);
    
    
    
    %   Delete NaN
    [I1,J] = find(isnan(xyzPoints));
    xyzPoints(I1,:) = [];
    
    %   Delete points with Z==0
    [I2,J] = find(xyzPoints(:,3)<1);
    xyzPoints(I2,:) = [];
    
    
    if(nargin > 5)
        xyzWeight = reshape(weight,[image_row*image_col,1]);
        xyzWeight(I1,:) = [];
        xyzWeight(I2,:) = [];
    end
%     %   For test
%     scatter3(xyzPoints(:,1),xyzPoints(:,2),xyzPoints(:,3),'.');
%     hold on;
%     quiver3(xyzPoints(:,1),xyzPoints(:,2),xyzPoints(:,3),xyzPoints(:,7),xyzPoints(:,8),xyzPoints(:,9));
end
