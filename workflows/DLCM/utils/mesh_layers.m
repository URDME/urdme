function layers = mesh_layers(N,nlayers)
%MESH_LAYERS Determine mesh boundary layers.
%   LAYERS = MESH_LAYERS(N,nlayers) given the neighbor matrix N,
%   returns nlayers of layers in the nlayers-by-2 cell-matrix
%   LAYERS. The first column contains descriptive names of the second
%   column, which contains indices into the mesh.
%
%   Example:
%     for mesh = 1:2
%       [P,E,T] = basic_mesh(mesh,20);
%       [V,R] = mesh2dual(P,E,T,'voronoi');
%       N = dt_neighe(V,R); N = (N ~= 0);
%       layers = mesh_layers(N,11);
%       cmap = colormap(hsv(size(layers,1)));
%         figure(mesh), clf,
%       for l = 1:size(layers,1)
%         patch('Faces',R(layers{l,2},:), ...
%               'Vertices',V, ...
%               'FaceColor',cmap(l,:), ...
%               'EdgeColor','black');
%         hold on,
%         axis([-1 1 -1 1]); axis equal
%       end
%     end
%
%   See also DT_OPERATORS, DT_NEIGHE, BASIC_MESH.

% S. Engblom 2024-10-28 (Simplified syntax) 
% S. Engblom 2024-09-30
  
% peel off the layers starting with corners and edges
layers = cell(nlayers,1);
neigh = N*ones(size(N,2),1);
% heuristics corner/"the middle"/edges:
ncorner = min(neigh);
nmiddle = max(neigh);
nedg = nmiddle-1;
% start the recursion:
layers{1} = sparse(neigh <= ncorner);
layers{2} = sparse(neigh <= nedg);
layers{2}(layers{1}) = 0;
visited = layers{1} | layers{2};
for i = 3:nlayers
  % neighbors of previous layer but not a member of *any* previous layers:
  layers{i} = N*layers{i-1} > 0;
  layers{i}(visited) = 0;
  visited = visited | layers{i};
end
% switch to index format
layers = cellfun(@find,layers,'UniformOutput',false);

% add labels as the first column
layers = [cellstr(reshape(sprintf('layer %2d',0:nlayers-1),[],nlayers)') ...
          layers];
layers{1} = 'corners';
layers{2} = 'edges';
