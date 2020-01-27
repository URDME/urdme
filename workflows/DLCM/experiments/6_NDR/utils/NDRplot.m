function NDRplot(D,U,V,R,mx)
%NDRplot Visualize Notch-Delta-Reporter pattern.
%   NDRplot(D,U,V,R,mx) visualizes the delta expression D in populated
%   voxels U over the (dual) mesh (V,R). The normalization used is to
%   scale D by mx, which is taken to be max(D(:)) if not given.
%
%   D is Nvoxels-by-M, where only the first column of D is
%   visualized. U is sparse Nvoxels-by-1 and counts the number of
%   cells in each voxel.
%
%   (V,R) is the mesh used for visualization, see MESH2DUAL.
%
%   No error-checking is performed.
%
%   See also BASIC_MESH, MESH2DUAL.

% S. Engblom 2018-01-24

% background mesh
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,

% normalization
if nargin < 5
  mx = full(max(D(:)));
end

% "white"
patch('Faces',R(U & D(:,1) <= 0.005*mx,:),'Vertices',V, ...
      'FaceColor',[254 254 254]/255,'EdgeColor',[0.7 0.7 0.7]);
% "orange"
patch('Faces',R(U & 0.005*mx < D(:,1) & ...
                D(:,1) <= 0.05*mx,:), ...
      'Vertices',V,'FaceColor',[252 108 15]/255,'EdgeColor',[0.7 0.7 0.7]);
% "light brown"
patch('Faces',R(U & 0.05*mx < D(:,1) & ...
                D(:,1) <= 0.5*mx,:), ...
      'Vertices',V,'FaceColor',[207 24 20]/255,'EdgeColor',[0.7 0.7 0.7]);
% "dark brown"
patch('Faces',R(U & 0.5*mx < D(:,1),:), ...
      'Vertices',V,'FaceColor',[62 11 1]/255,'EdgeColor',[0.7 0.7 0.7]);

axis equal, axis([-1 1 -1 1]);
set(gca,'xtick',[],'ytick',[]);
set(gca,'visible','off');
drawnow;
