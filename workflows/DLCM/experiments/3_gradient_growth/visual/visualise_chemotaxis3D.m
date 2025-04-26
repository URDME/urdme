%VISUALISE_CHEMOTAXIS Visualise chemotaxis model in 3D.

% E. Blom 2024-12-02

load chemotaxis_3D.mat

%% (1) Plot slices of chemotactic sensitivity
figure()
t = tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'tight');

idx = find(strcmp([umod.solverargs{:}], 'mumod')) + 1; % find mumod
[X, Y, Z] = meshgrid(-1:0.01:1,-1:0.01:1, -1:0.01:1);

str = ['a)', 'b)', 'c)', 'd)'];
n = 1;  % loop counter

for tt = [1 100]
  nexttile
  data = umod.solverargs{idx}{:}.U(1,:,tt);  % Chemical sensitivity
  UU = griddata(P(1,:), P(2,:), P(3,:), data, X, Y, Z);
  xslice = [-0.25,0,0.25]; yslice = []; zslice = [];
  s = slice(X,Y,Z,UU,xslice,yslice,zslice);
  s(1).EdgeColor = 'none';
  s(2).EdgeColor = 'none';
  s(3).EdgeColor = 'none';
  text(-0.3,+1,1.5, "$"+str(n:n+1)+"$", 'Interpreter','latex')
  n = n+2;
end

%set(gcf,'PaperPositionMode','auto');
%set(gcf,'Position',[100 100 340 160]);
%set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

% tweak these parameters to fully remove figure whitespace...
%t.Position = [0 -0.08 1 1.18]; % left, bottom, right, top?

% uncomment to save:
%exportgraphics(t,'chemtax3D_slices.pdf')

%% (2) plot cell distribution as points using standard color scheme

for tt = [1 100]
  nexttile
  adof = find(umod.U(1,:,tt)>0);
  %pdeplot3D(P, T, 'FaceAlpha', 0.2, ColorMapData=umod.U(1,:,1)) % domain
  srf = pdeplot3D(P, T, 'FaceAlpha', 0*0.1, 'EdgeColor', ...
    [0.3010 0.7450 0.9330]);
  srf.EdgeAlpha = 0.2;
  hold on
  plot3(P(1,adof),P(2,adof),P(3,adof), '.', 'color', ...  % singly occ.
    graphics_color('bluish green'),'markersize',14)
  sdof = find(umod.U(1,:,tt)>1);
  plot3(P(1,sdof),P(2,sdof),P(3,sdof), '.', 'color', ...  % doubly occ.
    graphics_color('vermillion'),'markersize',14)
  % cell 'shadows' (positions projected onto xy-plane)
  plot3(P(1,adof),P(2,adof), -ones(numel(adof),1), '.', ...
    'color', graphics_color('bluish green').*0.0,...
    'markersize',2)
  %plot3(P(1,sdof),P(2,sdof), -ones(numel(sdof),1), '.', ...
  %  'color', graphics_color('vermillion').*0.0, ...
  %  'markersize',2)
  %theta = -2*pi:0.05:2*pi;
  %plot3(cos(theta),sin(theta), -zeros(numel(theta),1), ...
  %  'k--', 'LineWidth',1.5)  % reference circles
  colorbar off
  view([45,-45, 20])
  axis(0.8.*[-1.0 1.0, -1.0, 1.0, -1.0, 1.0])
  lightangle(0,30)
  text(-0.5,-1.25,1.2, "$"+str(n:n+1)+"$", 'Interpreter','latex')
  n = n+2;
end

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 2*160]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

% tweak these parameters to fully remove figure whitespace...
t.Position = [0.052 0.0028 0.8983 0.9972]; % left, bottom, right, top?

% uncomment to save:
%exportgraphics(t,'chemtax3D_cells.pdf')