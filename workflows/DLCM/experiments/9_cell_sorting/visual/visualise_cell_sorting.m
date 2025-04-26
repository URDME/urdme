%VISUALISE_CELL_SORTING Visualise cell sorting dynamics.
%   Plots snapshots from simulation data

% E. Blom 2024-11-19

%% Exp 1, 1st row:
load cellsort1.mat

U = umod.U;

zoomf = 0.75; % domain zoom-in factor
nrows = 3;
str = ['a)', 'b)', 'c)', 'd)'];
n = 1;  % loop counter
figure()
t = tiledlayout(nrows,4, 'Padding', 'none', 'TileSpacing', 'none');
for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1+0.2 1+0.2])
  text(-0.7,+0.6,0, "$"+str(n:n+1)+"$", 'Interpreter','latex')
  axis off
  n = n+2;
end

%% Exp 2, 2nd row
load cellsort2.mat

U = umod.U;

for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1 1])
  plot(-1:1, 0.65*[1 1 1], 'k--', 'linewidth', 1.0)
  axis off
end

%% Exp 3, 3nd row
load cellsort3.mat

U = umod.U;

for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1-0.2 1-0.2])
  axis off
end

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 80*nrows]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

% tweak these parameters to fully remove figure whitespace...
t.Position = [0 -0.05 1 1.09]; % left, bottom, right, top?
% 2-by-4: t.Position = [0 -0.04 1 1.08];
% 3-by-4: t.Position = [0 -0.03 1 1.06];

% uncomment to save:
%exportgraphics(t,'cellsort1.pdf')