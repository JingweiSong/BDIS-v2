function xyzPoints = depth2XYZcamera(K, depth,image_row,image_col)

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
    
    %   Delete NaN
    [I,J] = find(isnan(xyzPoints));
    xyzPoints(I,:) = [];
end
