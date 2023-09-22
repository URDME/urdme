function [X,Y,Z] = meshdata(x,y,z)
% MESHDATA Converts 1-by-N^2 array of data, z, to N-by-N matrix, Z,
% to be used in MATLAB's meshgrid functions. The corresponding 1-by-N^2
% array positions x and y are converted in a compatible manner to X and Y
% such that the data Z(i,j) corresponds to the position [X(i,j), Y(i,j)].
%
%   input:  x - x-coordinate [1xN^2]
%           y - y-coordinate [1xN^2] 
%           z - data corresponding to (x,y) [1xN^2] 
%   output: X, Y, Z - NxN coordinates and data, corresponding to input
%
%   See also: MESHGRID

% E. Blom 2022-12-01

N = length(x);

if N ~= length(y)
    error('length of x not equal to length of y');
end
if N ~=  length(z)
    error('length of x not equal to length of data, z');
end

if height(x) ~= 1
    error('x dimensions not equal to (1,~)')
end

n = sqrt(N);

if round(n) - n ~= 0
    error('length of x is not a square')
end


X = zeros(n,n);
Y = zeros(n,n);
Z = zeros(n,n);

% re-index position arrays to matrices
for i = 1:n
    X(i,:) = x(1,((i-1)*n+1):i*n);
    Y(i,:) = y(1,((i-1)*n+1):i*n);
    Z(i,:) = z(1,((i-1)*n+1):i*n);
end

end