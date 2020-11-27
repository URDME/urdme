% Test M_Gamma function

[p,e,t] = simple_mesh(2);

[L,dM,N] = dt_operators(p,t);

% Mgamma-function
np = size(p,2); % Number of nodes 
nt = size(t,2); % Number of boundary edges
Mgamma = zeros(np,np); % Allocate mass matrix 

for N = 1:np
    
    % Find all edges that are connected to N
    
    % Find column number (triangle) in which the node N is part of 
    [~, col] = find(t(1:3,:)==N);
    
    % Find other nodes in these triangles
    tri = t(1:3, col);
    N_neig = unique(tri(tri~=N));
    N_rel = N_neig(N_neig > N); % Nodes with greater index and thus not accounted for
    
    if isempty(N_rel)
        E_rel = ones(1,3);
    else
        E_rel = ones(1,3).*sqrt((p(1,N_rel)-p(1,N)).^2 + (p(2,N_rel)-p(2,N).^2));
    end
    
    MK = [2 1 1; 1 2 1; 1 1 2].*E_rel/6; % element mass matrix
    Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb) + MK; % add element masses to M
    
end

Mgamma