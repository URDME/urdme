function [L,dM,N,Ax,Ay,Az] = dt_operators(P,T)
%DT_OPERATORS Operator assembly over (Delaunay) triangulation.
%   [L,dM,N] = DT_OPERATORS(P,T) assembles minus the Laplacian L over
%   the triangulation (P,T), along with the vector of voxel volumes
%   dM, and the sparse neighbor matrix N over the same mesh. The
%   function is built as an interface to PDE Toolbox assembly.
%
%   [L,dM,N,Ax,Ay,Az] = DT_OPERATORS(P,T) optionally assembles the
%   gradient/divergence operators Ax, Ay, Az.
%
%   Example uses:
%       %% (0) Create geometry get operators   
%       [P,~,T] = initmesh("lshapeg");
%       [L,dM,N,Ax,Ay,Az] = dt_operators(P,T);
%
%       %% (1) Solve Poisson equation
%       % define boundary
%       L(1,:) = 0;
%       L(1,1) = 1;               % only @ one node
%       rhs = ones(size(L,1),1);  % bc: u(1) = 1, sources = 1 everywhere
%       u = L\rhs;                % lumped mass matrix accounted for in L
%
%       %% (2a) get gradient of solution...
%       [L,M] = assema(P,T,1,1,0);
%       gradu(:,1) = M\(Ax*u);
%       gradu(:,2) = M\(Ay*u);    % might be non-smooth!
%       %% (2b) ...and divergence of the gradient:
%       divu(:,1) = M\(Ax*gradu(:,1)) + M\(Ay*gradu(:,2));
%
%   See also ASSEMA.

% E. Blom 2024-07-30 (addition, grad/div operators)
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

% optionally assemble advection matrices for projection of gradients
% and divergence (cf. Larson, M. G., and F. Bengzon. 2013. §10.1.3)
if nargout > 3
  np = size(P,2);
  nt = size(T,2);

  % essentially, component-wise convection matrix with advection velocity = 1
  Ax = sparse(np,np);    % x-component gradient operator
  Ay = sparse(np,np);    % y-component gradient operator
  Az = sparse(np,np);    % z-component gradient operator
  ndim = size(P,1);      % check whether to use 2D or 3D operators

  for i = 1:nt
    loc2glb = T(1:ndim+1,i);  % 2D - triangles, 3D - tetrahedra
    x = P(1,loc2glb);
    y = P(2,loc2glb);
    z = NaN;
    if ndim == 3
      z = P(3,loc2glb);
    end
    [area,b,c,d] = l_hat_grad(x,y,z);
    XK = ones(ndim+1,1).*b'.*area/(ndim+1);   % identical columns
    Ax(loc2glb,loc2glb) = Ax(loc2glb,loc2glb)+XK;
    YK = ones(ndim+1,1).*c'.*area/(ndim+1);   % identical columns
    Ay(loc2glb,loc2glb) = Ay(loc2glb,loc2glb)+YK;
    if ndim == 3
      ZK = ones(ndim+1,1).*d'.*area/(ndim+1); % identical columns
      Az(loc2glb,loc2glb)=Az(loc2glb,loc2glb)+ZK;
    end
  end
end
end
% ----------------------------------------------------------------------
function [area,b,c,d] = l_hat_grad(x,y,z)
%L_HAT_GRAD Constant gradients of the hat functions.
%   Adapted from M.G Larson, F. Bengzon, "The Finite Element Method:
%   Theory, Implementation, and Applications" (2013)

if isnan(z) % 2D
  area = polyarea(x,y);
  b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
  c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
  d = NaN;
else        % 3D
  [~,area] = convhull(x,y,z);   % volume
  b(1) = -(y(2)*(z(3)-z(4))+y(3)*(-z(2)+z(4))+y(4)*(z(2)-z(3)))/6/area;
  b(2) = -(y(1)*(z(3)-z(4))+y(3)*(-z(1)+z(4))+y(4)*(z(1)-z(3)))/6/area;
  b(3) = -(y(1)*(z(2)-z(4))+y(2)*(-z(1)+z(4))+y(4)*(z(1)-z(2)))/6/area;
  b(4) = -(y(1)*(z(2)-z(3))+y(2)*(-z(1)+z(3))+y(3)*(z(1)-z(2)))/6/area;

  c(1) = -(z(2)*(x(3)-x(4))+z(3)*(-x(2)+x(4))+z(4)*(x(2)-x(3)))/6/area;
  c(2) = -(z(1)*(x(3)-x(4))+z(3)*(-x(1)+x(4))+z(4)*(x(1)-x(3)))/6/area;
  c(3) = -(z(1)*(x(2)-x(4))+z(2)*(-x(1)+x(4))+z(4)*(x(1)-x(2)))/6/area;
  c(4) = -(z(1)*(x(2)-x(3))+z(2)*(-x(1)+x(3))+z(3)*(x(1)-x(2)))/6/area;

  d(1) = -(x(2)*(y(3)-y(4))+x(3)*(-y(2)+y(4))+x(4)*(y(2)-y(3)))/6/area;
  d(2) = -(x(1)*(y(3)-y(4))+x(3)*(-y(1)+y(4))+x(4)*(y(1)-y(3)))/6/area;
  d(3) = -(x(1)*(y(2)-y(4))+x(2)*(-y(1)+y(4))+x(4)*(y(1)-y(2)))/6/area;
  d(4) = -(x(1)*(y(2)-y(3))+x(2)*(-y(1)+y(3))+x(3)*(y(1)-y(2)))/6/area;
end
end
% ----------------------------------------------------------------------

