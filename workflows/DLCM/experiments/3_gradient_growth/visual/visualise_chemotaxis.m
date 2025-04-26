%VISUALISE_CHEMOTAXIS Visualise samples of various chemotaxis models.

% E. Blom 2024-12-03

%% (1) Plot slices of chemotactic sensitivity
figure()
t = tiledlayout(1,4, 'Padding', 'none', 'TileSpacing', 'none');
str = ['a)', 'b)', 'c)', 'd)'];
n = 1;  % loop counter
cmax = 2;

load chemotaxis_pressure.mat
nexttile
patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(U_samples), ...
'FaceColor','flat', 'EdgeColor','none');
clim([0,cmax])   % same color axis
hold on
theta = -2*pi:0.01:2*pi;
plot(0.25*cos(theta), 0.25*sin(theta), 'k', 'linewidth', 1)
axis([-0.98,0.3, -0.64, 0.64])
text(-0.9,+0.5,0, "$"+str(n:n+1)+"$", 'Interpreter','latex')
n = n+2;
axis square
axis off

load chemotaxis_diffusion.mat
nexttile
patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(U_samples), ...
'FaceColor','flat', 'EdgeColor','none');
clim([0,cmax])   % same color axis
hold on
plot(0.25*cos(theta), 0.25*sin(theta), 'k', 'linewidth', 1)
axis([-0.98,0.3, -0.64, 0.64])
text(-0.9,+0.5,0, "$"+str(n:n+1)+"$", 'Interpreter','latex')
n = n+2;
axis square
axis off

%load chemotaxis_cons3_attr.mat
load chemotaxis_cons.mat
nexttile
patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(U_samples), ...
'FaceColor','flat', 'EdgeColor','none');
clim([0,cmax])   % same color axis
hold on
plot(0.25*cos(theta), 0.25*sin(theta), 'k', 'linewidth', 1)
axis([-0.98,0.3, -0.64, 0.64])
text(-0.9,+0.5,0, "$"+str(n:n+1)+"$", 'Interpreter','latex')
n = n+2;
axis square
axis off

load chemotaxis_modLa.mat
nexttile
patch('Faces',R(:,:),'Vertices',V,'FaceVertexCData',full(U_samples), ...
'FaceColor','flat', 'EdgeColor','none');
clim([0,cmax])   % same color axis
hold on
plot(0.25*cos(theta), 0.25*sin(theta), 'k', 'linewidth', 1)
axis([-0.98,0.3, -0.64, 0.64])
text(-0.9,+0.5,0, "$"+str(n:n+1)+"$", 'Interpreter','latex')
n = n+2;
axis square
axis off
colorbar

% use the traditional dlcm cell colors in continuous format
mg = [0.9 0.9 0.9];     % i) 'gray'
mb = [0 158 115]./255;  % ii) 'bluish green'
mr = [213 94 0]./255;   % iii) 'vermillion' -- cf. graphics_color
map = zeros(256,3);     % colormap i)-ii)-iii):
for i = 1:3
  map(1:128,i) = linspace(mg(i),mb(i), 128);
end
for i = 1:3
  map(129:256,i) = linspace(mb(i),mr(i), 128);
end
colormap(map)

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 90]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

% tweak these parameters to fully remove figure whitespace...
t.Position = [0 -0.15 0.885 1.32]; % left, bottom, right, top?

% uncomment to save:
%exportgraphics(t,'chemtax_variety.pdf')