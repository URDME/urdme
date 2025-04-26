%VISUALISE_DELTA_NOTCH Visualise Delta-Notch dynamics
%   Plots snapshots of the dynamics for both continuous and discrete
%   models

% E. Blom 2024-11-08

%% Continuous DN dynamics
load contDNgrowth.mat

idx = find(strcmp([umod.solverargs(:)], 'mumod')) + 1; % find mumod
X = reshape(umod.solverargs{idx}.U(1,:,:), ...
size(umod.U,2), size(umod.U,3));             % = mumod.u0(1,:,:)

Smax = 400;   % maximum value
X = X./Smax;  % normalize to [0,1]

figure()
t = tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'none');
tspan = [1 10 25 40];

n = 1;      % tile counter
for i = tspan
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(X(:,i)>0);
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  clim([0 1])
    xmax = 1;
  if n == 4
    xmax = 1;   % extend final tile per row
  end
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

%% Discrete DN dynamics
load discDNgrowth.mat

state = 3;                             % which state to vis.

idx = find(strcmp([umod.solverargs(:)], 'mumod')) + 1; % find mumod
X = reshape(umod.solverargs{idx}.U(1,:,:), ...
size(umod.U,2), size(umod.U,3));             % = mumod.u0(1,:,:)

Smax = 400;   % maximum value
X = X./Smax;  % normalize to [0,1]

n = 1;        % tile counter
for i = tspan
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(X(:,i)>0);
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  clim([0 1])
  xmax = 1;
  if n == 4
    xmax = 1;   % extend final tile per row
  end
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 360 190]);
set(gca, 'fontname', 'Roman', 'FontSize', 8.0)

% create and use colormap akin to the paraview standard
m1 = [0.1392 0.2083 0.6469];    % i) blue
m2 = [0.5096 0.0458 0.0406];    % ii) red
mw = [0.8 0.8 0.8];             % iii) white-grey
map = zeros(256, 3);            % create colormap i)-iii)-ii):
for i = 1:3
  map(1:128,i) = linspace(m1(i),mw(i), 128);
end
for i = 1:3
  map(128:256,i) = linspace(mw(i),m2(i), 129);
  map(128:256,3) = map(128:256,3) - linspace(0,0.03,129)'; % adjust
  map(128:256,2) = map(128:256,2) + linspace(0,0.05,129)';
end
map = map + 0.2*[0.9 0.9 0.9];  % brighten
colormap(map)

cb = colorbar;
set(cb,'Position',[0.91 0.05 .015 0.9])

% tweak these parameters to fully remove figure whitespace...
t.Position = [-0.09 0.00 1 1.00]; % left, bottom, right, top

% uncomment to save:
%exportgraphics(t,'DN.pdf', 'resolution', 600)