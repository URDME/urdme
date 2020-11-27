% Function that allocates the boundary matrix Mgamma.
% Output: M the allocated mass matrix
% Input: p the point matrix, t the connectivity matrix 
function Mgamma = assemble_Mgamma2(p,t,imp)
    np = size(p,2); % Number of nodes 
    nt = size(t,2); % Number of boundary edges
    Mgamma = sparse(np,np); % Allocate mass matrix 
    inds = [1,2;2,3;3,1];
if imp == 1
    tic
    for K = 1:nt
        loc2glb = t(1:3,K); 
        x = p(1,loc2glb); % x-coordinates of triangle nodes
        y = p(2,loc2glb); % y-coordinates of triangle nodes
        for ii = 1:length(loc2glb)
           len = sqrt((x(inds(ii,1))-x(inds(ii,2)))^2 ...
               +(y(inds(ii,1))-y(inds(ii,2)))^2);
           MK = [2 1; 1 2]/(2*6)*len; % element mass matrix
           Mgamma(loc2glb(inds(ii,:)),loc2glb(inds(ii,:))) = ...
               Mgamma(loc2glb(inds(ii,:)),loc2glb(inds(ii,:)))+ MK;
        end
    end
    toc
elseif imp == 2
    % Other implementations
    len = zeros(1,3);
    tic
    for K = 1:nt
        loc2glb = t(1:3,K); 
        x = p(1,loc2glb); % x-coordinates of triangle nodes
        y = p(2,loc2glb); % y-coordinates of triangle nodes
        for ii = 1:length(loc2glb)
           len(1,ii) = sqrt((x(inds(ii,1))-x(inds(ii,2)))^2 ...
               +(y(inds(ii,1))-y(inds(ii,2)))^2);
        end
        MK = ([2 0 1; 1 2 0; 0 1 2].*len + ...
            [2 1 0; 0 2 1; 1 0 2].*circshift(len,1))/(2*6); % element mass matrix 
        Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb)+ MK;
    end
    toc
else
    disp('No such Mgamma assembly implementation, returning empty matrix');
end
end