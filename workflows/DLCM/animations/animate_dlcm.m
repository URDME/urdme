%ANIMATE_DLCM DLCM animation script examples.
%   How to animate the different variables of a DLCM-simulation using
%   the umod-struct. Simply copy and paste a section.

% E. Blom 2024-12-01

% Dual mesh needed, e.g., [V,R] = mesh2dual(P,E,T,'voronoi');

%% (1) animating number of cells:
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
for i = 1:numel(umod.tspan)
  patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(U(:,i)), ...
  'FaceColor','flat', 'EdgeColor','none');
  drawnow
end

%% (2) animating the cell types in one layer of cells:
X = umod.U(:,:,:);
figure();
for i = 1:numel(umod.tspan)
   clf,
  for n = 1:ntypes
    hold on
    patch('Faces',R(X(n,:,i)>0,:),'Vertices',V, ...
     'FaceColor',0.5/n*[1,1,1], 'EdgeColor','none');
  end
  drawnow
end

%% (3) animate cell internal states
idx = find(strcmp([umod.solverargs(:)], 'mumod')) + 1; % find mumod
data = reshape(umod.solverargs{idx}.U(1,:,:), ...
  size(umod.U,2), size(umod.U,3));             % = mumod.u0(1,:,:)
for i = 1:numel(umod.tspan)
  patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(data(:,i)), ...
  'FaceColor','flat', 'EdgeColor','none');
  drawnow
end


%% (4) animate micro-environment quantities
quant = 1;        % which quantity to plot
Q = reshape(umod.U(ntypes + quant,:), ...
size(umod.U,2),size(umod.U,3));
for i = 1:numel(umod.tspan)
  patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(Q(:,i)), ...
    'FaceColor','flat', 'EdgeColor','none');
  drawnow
end

%% (5) plot cell trajectories
Nstates = size(umod.U,1)/2;
glix = umod.private.dlcm.glix;
for t = 2:numel(umod.tspan)      % map indices to voxels versus time
  for k = 2:max(glix(2,:,1))     % find all cells up to this index
    idx(t,k) = find(glix(1,:,t)==k | glix(2,:,t)==k); 
  end
end
figure(1), clf, hold on,
K = min(max(glix(2,:,1)),30);
for k = 2:K % plot a few of cell trajectories
plot(P(1, idx(2:end,k)), P(2,idx(2:end,k)), 'Color', k/K*[1,1,1]);
plot(P(1, idx(end,k)), P(2,idx(end,k)), '*', 'Color', k/K*[1,1,1]);
end
axis([-1, 1, -1, 1])

%% (6) 3D visualisation

% most 3D vis. tools require the 'meshgrid'-format:
[X, Y, Z] = meshgrid(-1:0.01:1,-1:0.01:1, -1:0.01:1);
t = 100;     % tt = umod.tspan(t)
% interpolate data to this X,Y,Z grid (e.g., pressure, cell number, cell
% states...)
idx = find(strcmp([umod.solverargs(:)], 'mumod')) + 1; % find mumod
state = 1;
data = umod.solverargs{idx}.U(state,:,t);   % internal state 1
% data = umod.U(2,:,t)  % Pressure
UU = griddata(P(1,:), P(2,:), P(3,:), data, X, Y, Z);

% now use whichever tools you want:

% (6a) slicing
figure()
xslice = [-0.25,0,0.25]; yslice = []; zslice = [];
s = slice(X,Y,Z,UU,xslice,yslice,zslice);
s(1).EdgeColor = 'none';
%s(1).FaceAlpha = 0.75;
s(2).EdgeColor = 'none';
%s(2).FaceAlpha = 0.75;
s(3).EdgeColor = 'none';
%s(3).FaceAlpha = 0.75;

% (6b) contours of above slices + contours in xy-plane:
figure()
contourslice(X,Y,Z,UU,xslice,yslice,zslice);
contourslice(X,Y,Z,UU,[],[],0);
view(3)

% (6c) isosurface
figure()
isosurface(X,Y,Z,UU,0.04)

% (6d) plot cell distribution as points using standard color scheme
figure()
adof = find(umod.U(1,:,t)>0);
pdeplot3D(P, T, 'FaceAlpha', 0.2, ColorMapData=umod.U(1,:,1)) % domain
hold on
plot3(P(1,adof),P(2,adof),P(3,adof), '.', 'color', ...  % singly occ.
  graphics_color('bluish green'),'markersize',10)
sdof = find(umod.U(1,:,t)>1);
plot3(P(1,sdof),P(2,sdof),P(3,sdof), '.', 'color', ...  % doubly occ.
  graphics_color('vermillion'),'markersize',10)
% cell 'shadows' (positions projected onto xy-plane)
plot3(P(1,adof),P(2,adof), -ones(numel(adof),1), 'k.', 'markersize',7)
