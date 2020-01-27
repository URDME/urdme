function [P,E,T,gradquotient,V,R,L,dM,N,Nj,Np,U,param] = ...
    DLCMlayer(Nvoxels,Theta,dTheta,len)
%DLCMlayer Utility function which sets up the DLCM-layer.

% S. Engblom 2018-02-08 (minor revision)
% S. Engblom 2018-01-24 (minor revision)
% S. Engblom 2017-09-14

% repeatable
rng(456);

% build the hex-mesh
[P,E,T,gradquotient] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');
Nvoxels = size(P,2);

% radius taken to be max distance from all vertices to voxel center
Radius = R;
for j = 1:size(R,1)
  % all NaN's mapped to first vertex
  Radius(j,isnan(Radius(j,:))) = R(j,1);
end
Radius = sqrt(max(sum(tsum(reshape(V(Radius,:),Nvoxels,size(R,2),2),-P, ...
                           [1 3 2],[2 1]).^2,2),[],3));

% additional small stochastic perturbation of radii
%Radius = mean(Radius)*(1+0.05*randn(size(Radius)));

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);
Nj = N;
% N is neighbor matrix, Nj is the junctional contact neighbor matrix

% population of cells (constant here, but format allows for a growing
% domain)
U = sparse(ones(Nvoxels,1));

% *** known issue: Theta is taken to be static here

% protrusions in direction Theta with angle width dTheta and length
% len * cell-radius
if nargin < 2
  warning('Non-existing variable (Theta). Default used.');
  Theta = pi*sprand(U);
else
  if Theta > 0
    Theta = Theta*spones(U);
  elseif Theta == -1
    % special syntax: understands perpendicular to radii
    warning('Understands Theta perpendicular to radii.');
    Theta = sparse(atan2(P(1,:),-P(2,:))');
    % visualize:
% $$$       ang = atan2(P(1,:),-P(2,:));
% $$$       quiver(P(1,:),P(2,:),cos(ang),sin(ang)), axis equal
  elseif Theta == -2
    % special syntax: understands radial direction
    warning('Understands Theta in the radial direction.');
    Theta = sparse(atan2(P(2,:),P(1,:))');
    % visualize:
% $$$       ang = atan2(P(2,:),P(1,:));
% $$$       quiver(P(1,:),P(2,:),cos(ang),sin(ang)), axis equal
  elseif Theta == -3
    % special syntax: "spiral"
    Theta = sparse(pi/6+atan2(P(1,:),-P(2,:))');
    % visualize:
% $$$       ang = pi/6+atan2(P(1,:),-P(2,:));
% $$$       quiver(P(1,:),P(2,:),cos(ang),sin(ang)), axis equal
  else
    error('Could not parse Theta value.');
  end
end
if nargin < 3
  warning('Non-existing variable. Default used.');
  dTheta = pi;
end
if nargin < 4
  warning('Non-existing variable. Default used.');
  len = 3.5;
end
% (radius ~ 2 a.u., protrusion length = len*radius)

% static protrusional maximum neighbor matrix
if len ~= 0
  Npmax = connect(complex(P(1,:),P(2,:)),2*len*max(Radius));
else
  Npmax = sparse(size(Nj,1),size(Nj,2));
end

% build the "dynamic" (it is static in this experiment) neighbor
% matrix - for protrusional contacts
iip = zeros(0,1); jjp = zeros(0,1);
iij = zeros(0,1); jjj = zeros(0,1);
for i = find(U(:,1))'
  for j = find(Npmax(:,i))'
    if j ~= i
      % how is cell i affected by cell j?
      p = contact(P(:,i),Radius(i),len*Radius(i),full(Theta(i,1)), ...
                  P(:,j),Radius(j),len*Radius(j),full(Theta(j,1)), ...
                  dTheta);
      if p == 1
        % j --> i is protrusional: either the protrusions touch symmetrically
        % or the protrusion of j touches the membrane of i
        iip = [iip; i]; jjp = [jjp; j];
      elseif p == 2
        % j --> i is junctional: either the membranes of i and j touches
        % symmetrically or the membrane of j touches the protrusion
        % of i
        iij = [iij; i]; jjj = [jjj; j];         
      elseif p == 3
        % j --> i is both protrusional and junctional
        iip = [iip; i]; jjp = [jjp; j];
        iij = [iij; i]; jjj = [jjj; j];
      end
    else
      % i == j, the convention here is that cells in the same voxel implies
      % both a protrusional and a junctional contact
      iip = [iip; i]; jjp = [jjp; j];
      iij = [iij; i]; jjj = [jjj; j];
    end
  end
end
Njx = fsparse(iij,jjj,1,size(Npmax));
Nj = Nj+Njx;
Nj = spreplace(Nj,ones(nnz(Nj),1)); % (handle double contacts)
Np = fsparse(iip,jjp,1,size(Npmax));

% delta-notch-reporter model parameters, baseline:
param.betaN = 100;
param.betaD = 500;
param.betaR = 300000;

param.kt = 2;
param.kc = 0.5;
param.kRS = 1e7;

param.wa = 1;
param.qa = 1;
param.wb = 1;
param.qb = 1;
