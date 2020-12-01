function [P,E,T,gradquotient] = basic_mesh2(mesh,Nvoxels)
% BASIC_MESH Basic regular mesh for population-level cell models.
%   [P,E,T,gradquotient] = BASIC_MESH(mesh,Nvoxels) creates a basic
%   rectangular geometry discretized either in a Cartesian mesh (mesh
%   = 1) or in a hexagonal mesh (mesh = 2). The discretization is
%   Nvoxels-by-Nvoxels and is either square (mesh = 1) or rectangular
%   (mesh = 2).
%
%   On return, (P,E,T) is a PDE Toolbox triangulation and the scalar
%   gradquotient is the edge length to voxel center distance quotient.
%
%   Examples:
%     for mesh = 1:2
%       [P,E,T] = basic_mesh(mesh,10);
%       [V,R] = mesh2dual(P,E,T,'voronoi');
%       figure, patch('Faces',R,'Vertices',V, ...
%          'FaceColor',[0.9 0.9 0.9],'EdgeColor','black');
%       hold on,
%
%       % these are the voxel center points:
%       plot(P(1,:),P(2,:),'ko');
%
%       % this is the computational mesh:
%       triplot(T(1:3,:)',P(1,:),P(2,:),'b--');
%       axis([-1 1 -1 1]); axis equal
%     end
%
%   See also MESH2DUAL, PATCH, TRIPLOT.

% S. Engblom 2017-08-29

if mesh == 1
  % square Cartesian mesh
  gd = [3 4 -1 1 1 -1 -1 -1 1 1]';
  sf = 'SQ1';
  ns = char(sf)';
  G = decsg(gd,sf,ns);
  [P,E,T] = poimesh2(G,Nvoxels-1);

  % integration of pressure gradients over edges

  % distance between two voxel midpoints
  h = 2/(Nvoxels-1);
  D1 = h;

  % edge length
  L1 = h;

  % quotient
  gradquotient = L1/D1; % (= 1 for Cartesian mesh)

elseif mesh == 2
  % hexagonal mesh
  
  % height and width
  h = 2/((Nvoxels-1)*sqrt(3)/2); % such that numel(xlayer) = Nvoxels
  w = sqrt(3)/2*h;

  % layers
  xlayer = (-1:w:1)'; % [-1,1]
  xlayer = xlayer-mean(xlayer)-w/4;
    
  lim = 3/4*2/sqrt(3); % [-lim,lim], such that numel(ylayer) = Nvoxels
  ylayer = -lim:3/4*h:lim;

  % mesh points
  x = repmat(xlayer,[1 numel(xlayer)])+ ...
      repmat(w/2*mod((0:numel(ylayer)-1),2),[numel(xlayer) 1]);
  y = repmat(ylayer,[numel(xlayer) 1]);
  P = [x(:)'; y(:)'];

  % triangulate the points
  DT = delaunayTriangulation(P');
  T = DT.ConnectivityList';
  T = [T; ones(1,size(T,2))];

  % build edge matrix
  E = convexHull(DT);
  E = [E(1:end-1)'; E(2:end)'];
  E = [E; zeros(1,size(E,2)); ones(2,size(E,2)); ...
       ones(1,size(E,2)); zeros(1,size(E,2))];

  % integration of pressure gradients over edges

  % distance between two voxel midpoints
  D1 = sqrt((3/4*h)^2+(w/2)^2);

  % edge length
  L1 = sqrt((1/4*h)^2+(w/2)^2);

  % quotient
  gradquotient = L1/D1;
else
  error('Unsupported mesh case.');
end
