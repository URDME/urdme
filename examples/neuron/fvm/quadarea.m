% POLYAREA3D 
%
% Computes the normal and area of a 3D quadrilateral. 
% 

function [n,ar] = quadarea(X)

% diagonals of the quadrilateral
v1 = X(:,3)-X(:,1);
v2 = X(:,4)-X(:,2);

% Normal to the surface
n  = [v1(2)*v2(3)-v1(3)*v2(2),v1(3)*v2(1)-v1(1)*v2(3), ... 
      v1(1)*v2(2)-v1(2)*v2(1)];
  
ar = norm(n);
n  = n/ar;
% area
ar = 0.5*ar;
