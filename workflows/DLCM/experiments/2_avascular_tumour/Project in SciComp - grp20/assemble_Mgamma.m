function Mgamma = assemble_Mgamma(p,t)
%ASSEMBLE_MGAMMA allocates the boundary matrix Mgamma
%   Input: p the point matrix, t the connectivity matrix
%   Output: Mgamma the allocated mass matrix
%
%   The structure of each element MK is as follows:
%     [2*(len(1)+len(3)),     len(1),         len(3);
%      len(1),            2*(len(1)+len(2)),  len(2);
%      len(3),            len(2),         2*(len(2)+len(3))]


    np = size(p,2); % Number of nodes 
    nt = size(t,2); % Number of boundary edges
    Mgamma = sparse(np,np); % Allocate mass matrix 
    inds = [1,2;2,3;3,1]; % Indices for the points in each triangle
    len = zeros(1,3); % Lenght vector
    
    % Loop over all triangles
    for K = 1:nt
        loc2glb = t(1:3,K); % The nodes that spans triangle
        x = p(1,loc2glb); % x-coordinates of triangle nodes
        y = p(2,loc2glb); % y-coordinates of triangle nodes
        
        % Get all boundary lengths and save them in len
        for ii = 1:length(loc2glb)
           len(1,ii) = sqrt((x(inds(ii,1))-x(inds(ii,2)))^2 ...
               +(y(inds(ii,1))-y(inds(ii,2)))^2);
        end
        
        % Create MK element
        MK = ([2 0 1; 1 2 0; 0 1 2].*len + ...
            [2 1 0; 0 2 1; 1 0 2].*circshift(len,1))/(2*6);     
        % Add to Mgamma
        Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb)+ MK;
    end
end