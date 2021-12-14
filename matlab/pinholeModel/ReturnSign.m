function [ bool ] = ReturnSign( x1,y1,x2,y2,x3,y3,x,y )
%RETURNSIGN Summary of this function goes here
%   Detailed explanation goes here
%   x1,y1,x2,y2,x3,y3 are the coordinate of the triangle
%   x,y is the target point

% Judge whether the point falls in the triangle by the equation
% F(x,y) = y - A*x - B;
% A = (y2-y1)/(x2-x1) B =(x2*y1-x1*y2)/(x2-x1)
% A triangle C(x1,y1) D(x2,y2) E(x3,y3), If point P IS inside, C and P are of same sign for another line,
% So does D and P, E and P.

% For point D and target point
if(x3 ~= x2)
    A = (y3 - y2) / (x3 - x2);
    B = (x3 * y2 - x2 * y3) / ( x3 - x2 );
    F_D = y1 - A*x1 - B;
    F_D_P = y - A*x - B;
else  
    if((x1-x2)*(x-x2)>=0)
        F_D = 1;
        F_D_P = 1;
    else
        F_D = 1;
        F_D_P = -1;
    end
end

% For point E and target point
if(x3 ~= x1)
    A = (y3 - y1) / (x3 - x1);
    B = (x3 * y1 - x1 * y3) / ( x3 - x1 );
    F_E = y2 - A*x2 - B;
    F_E_P = y - A*x - B;
else
    if((x2-x1)*(x-x1)>=0)
        F_E = 1;
        F_E_P = 1;
    else
        F_E = 1;
        F_E_P = -1;
    end
end

% For point F and target point
if(x2 ~= x1)
    A = (y2 - y1) / (x2 - x1);
    B = (x2 * y1 - x1 * y2) / ( x2 - x1 );
    F_F = y3 - A*x3 - B;
    F_F_P = y - A*x - B;
else
    if((x3-x1)*(x-x1)>=0)
        F_F = 1;
        F_F_P = 1;
    else
        F_F = 1;
        F_F_P = -1;
    end
end

if((F_D * F_D_P >= 0) && (F_E * F_E_P >= 0) && (F_F * F_F_P >= 0))
    bool = 1;
else
    bool = 0;
end
end

