
% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[p,e,t,gradquotient] = basic_mesh(1,Nvoxels);

np = size(p,2); % Number of nodes 
nt = size(t,2); % Number of boundary edges
Mgamma = sparse(np,np); % Allocate mass matrix 

len = zeros(1,3);
for K = 1:nt
    loc2glb = t(1:3,K); % triangle nodes in triangle K
    x = p(1,loc2glb); % x-coordinates of triangle nodes
    y = p(2,loc2glb); % y-coordinates of triangle nodes
    for ii = 1:length(loc2glb)
       len(1,ii) = sqrt((x(1+mod(ii-1,length(loc2glb)))-x(1+mod(ii,length(loc2glb))))^2 ...
           +(y(1+mod(ii-1,length(loc2glb)))-y(1+mod(ii,length(loc2glb))))^2); 
    end
    MK = [2 1 1; 1 2 1; 1 1 2]/12.*len; % element mass matrix 
    Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb)+ MK; % add element masses to M
end
