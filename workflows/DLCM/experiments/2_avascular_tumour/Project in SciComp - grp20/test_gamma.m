% ---------------------------------------------------------------------- %
% Test down-scaled problem
% ---------------------------------------------------------------------- %

% Parameters 
Nvoxels = 2;
alpha = 10.^(-3:5);
alpha_inv = 1./alpha;

% Create mesh and assemble matrices
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% Multiply L by inverse of M
[L0,M] = assema(P,T,1,1,0);

% Invert L using lumped mass matrix
% the (lumped) mass matrix dM gives the element volume
dM = full(sum(M,2));
ndofs = size(dM,1);

% explicitly invert the lumped mass matrix and filter the diffusion matrix
[i,j,s] = find(L0);
s = s./dM(i);
keep = find(i ~= j); % (removes only the diagonal)
i = i(keep); j = j(keep); s = s(keep);

% rebuild L, ensuring that the diagonal equals minus the sum of the
% off-diagonal elements
LL = sparse(i,j,s,ndofs,ndofs);
LL = LL+sparse(1:ndofs,1:ndofs,-full(sum(LL,2)));

N = sparse(i,j,1,ndofs,ndofs);

% Assemble boundary matrix and filter it as the
Mgamma = assemble_Mgamma(P,T);
Mgamma_dM = Mgamma./dM;

% Create DOFs
adof = 1:floor(Nvoxels^2/2); % Active voxels (source voxels)
idof = adof(end)+1:size(LL,1); 

% Laplacian adjustment
Lai = fsparse(idof,idof,1,size(LL));

% Investigate matrix multiplication for different values of alpha
RHS = ones(size(LL,1),1) + full(fsparse(adof',1,1./dM(adof),[size(LL,1) 1]));
for a_inv = 1e6%alpha_inv
    
    % Case one: All matrices being inverted at the same time
    LHS1 = (M \ (L0 - Lai*L0 + a_inv*Lai*Mgamma));
    X1 = LHS1 \ RHS;
    
    % Case two: L and M inverted in dt_operators.m
    LHS2 = LL - Lai*LL + a_inv*Lai*Mgamma_dM;
    X2 = LHS2 \ RHS;
    
end

% ---------------------------------------------------------------------- %
% Test M_Gamma function
% ---------------------------------------------------------------------- %

% [p,e,t] = simple_mesh(2);
% 
% [L,dM,N] = dt_operators(p,t);
% 
% % Mgamma-function
% np = size(p,2); % Number of nodes 
% nt = size(t,2); % Number of boundary edges
% Mgamma = zeros(np,np); % Allocate mass matrix 

% for N = 1:np
%     
%     % Find all edges that are connected to N
%     
%     % Find column number (triangle) in which the node N is part of 
%     [~, col] = find(t(1:3,:)==N);
%     
%     % Find other nodes in these triangles
%     tri = t(1:3, col);
%     N_neig = unique(tri(tri~=N));
%     N_rel = N_neig(N_neig > N); % Nodes with greater index and thus not accounted for
%     
%     if isempty(N_rel)
%         E_rel = ones(1,3);
%     else
%         E_rel = ones(1,3).*sqrt((p(1,N_rel)-p(1,N)).^2 + (p(2,N_rel)-p(2,N).^2));
%     end
%     
%     MK = [2 1 1; 1 2 1; 1 1 2].*E_rel/6; % element mass matrix
%     Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb) + MK; % add element masses to M
%     
% end

% Mgamma