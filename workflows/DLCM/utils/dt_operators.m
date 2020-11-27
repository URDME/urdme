function [L,dM,N,M] = dt_operators(P,T)
%DT_OPERATORS Operator assembly over (Delaunay) triangulation.
%   [L,dM,N] = DT_OPERATORS(P,T) assembles minus the Laplacian L over
%   the triangulation (P,T), along with the vector of voxel volumes
%   dM, and the sparse neighbor matrix N over the same mesh. The
%   funtion is built as an interface to PDE Toolbox assembly.
%
%   See also ASSEMA.

% S. Engblom 2017-12-20 (revision, changed filtering)
% S. Engblom 2017-08-29

% assemble minus the Laplacian on this grid (ignoring BCs) as well as
% the Mass-matrix
[L,M] = assema(P,T,1,1,0);

% the (lumped) mass matrix gives the element volume
dM = full(sum(M,2));
ndofs = size(dM,1);

% explicitly invert the lumped mass matrix and filter the diffusion matrix
[i,j,s] = find(L);
s = s./dM(i);
%keep = find(s < 0); % (possibly removes negative off-diagonal elements)
keep = find(i ~= j); % (removes only the diagonal)
i = i(keep); j = j(keep); s = s(keep);

% rebuild L, ensuring that the diagonal equals minus the sum of the
% off-diagonal elements
L = sparse(i,j,s,ndofs,ndofs);
L = L+sparse(1:ndofs,1:ndofs,-full(sum(L,2)));

% static neighbor matrix
N = sparse(i,j,1,ndofs,ndofs);
