function [ depthimage ] = ProjectModel( face, vertices, num_imagerow,num_imagecolunm, K )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function: Points to depth(pinhole projection) 
%   Method:   
%   Input:    
%             face:             Faces of model
%             vertices:         Triangle mesh
%             num_imagerow:     Row of image
%             num_imagecol:     Column of image
%             K:                3¡Á3 intrinsic matrix
%   Returns:   
%             depthimage:       Depth image
%   Author:   Jingwei Song.     25/07/2016
%   PS:
%   face: Stores the index of each points constructing a face
%         1 2 3; 4 5 6; 7 8 9;
%   vertices: stores the position of each point
%         0 0 0; 4 4 4; 2 2 2;
%   Intersection between ray and a plane. See picture
%   "formulation_ray_plane"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

depthimage = zeros(num_imagerow,num_imagecolunm);
num_face = size(face,1);
for i = 1 : num_face
    %get three point on the face
    index_pointA = face(i,1);
    index_pointB = face(i,2);
    index_pointC = face(i,3);
    pointA = vertices(index_pointA,:);
    pointB = vertices(index_pointB,:);
    pointC = vertices(index_pointC,:);
    
    if((pointA(1) - pointB(1) == 0) && (pointA(2) - pointB(2) == 0))
        continue;
    end
    if((pointA(1) - pointC(1) == 0) && (pointA(2) - pointC(2) == 0))
        continue;
    end
    if((pointC(1) - pointB(1) == 0) && (pointC(2) - pointB(2) == 0))
        continue;
    end
    
    triangle = zeros(3,3);      %   column1:u column2:v column3:depth value
    [triangle(:,1),triangle(:,2),triangle(:,3)] = ...
                XYZcamera2depth(K,[pointA;pointB;pointC]);
            
    max_x = max(triangle(:,1));
    max_y = max(triangle(:,2));
    max_z = max(triangle(:,3));
    min_x = min(triangle(:,1));
    min_y = min(triangle(:,2));
    min_z = min(triangle(:,3));
    
    
    %   If all of the triangle is outside of the image, ignore the face
    if ((max_x<=0)||(min_x>num_imagerow)||(max_y<=0)||(min_y>num_imagecolunm))
        continue;
    end
    if(max_x>num_imagerow)
        max_x = num_imagerow;
    end
    if(min_x<1)
        min_x = 1;
    end
    if(max_y>num_imagecolunm)
        max_y = num_imagecolunm;
    end
    if(min_y<1)
        min_y = 1;
    end
%     if ((min_x<=0)||(max_x>num_imagerow)||(min_y<=0)||(max_y>num_imagecolunm))
%         continue;
%     end
    
    % triangular function ax+by+cz+d=0
    % z=(-d-ax-by)/c   Calculate Z value by the surface equation
    x1=pointA(1);x2=pointB(1);x3=pointC(1);y1=pointA(2);y2=pointB(2);y3=pointC(2);z1=pointA(3);z2=pointB(3);z3=pointC(3);
    a = (y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1);
    b = (z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1);
    c = (x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1);
    d = - (a * x1 + b * y1 + c * z1);
    for m = min_x : max_x
        for n = min_y : max_y
            if(c == 0)
                continue;
            end
            z_point = -d/(...
                          a*(n-K(1,3))/K(1,1)+b*(m-K(2,3))/K(2,2)+c...
                         );
            x_point = (n-K(1,3))/K(1,1)*z_point;
            y_point = (m-K(2,3))/K(2,2)*z_point;
            bool = ReturnSign( pointA(1),pointA(2),pointB(1),pointB(2),pointC(1),pointC(2),x_point,y_point);
            if(bool > 0)
                                
                if(depthimage(m,n) == 0)
                    depthimage(m,n) = z_point;
                else
                    if(z_point < depthimage(m,n))
                        depthimage(m,n) = z_point;                
                    end
                end
            end
        end
    end
   
    
end

end

