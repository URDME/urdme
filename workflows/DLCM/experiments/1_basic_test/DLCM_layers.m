%DLCM_LAYERS Visualizes DLCM multi-level structure.

% E. Blom 2024-12-17

% get hexagonal mesh
N = 8;     % N^2 = number of voxels
[P,E,T,~] = basic_mesh(2,N);
[V,R] = mesh2dual(P,E,T,'voronoi'); 

% plot big hexagon...
figure(1), clf,
t = tiledlayout(1, 1, 'tilespacing', 'none', 'padding', 'tight');
nexttile
h = 0.7;
plot3([-sqrt(3/4)*h 0 sqrt(3/4)*h sqrt(3/4)*h 0 -sqrt(3/4)*h ...
  -sqrt(3/4)*h], [h/2 h h/2 -h/2 -h -h/2 h/2], -0.5*ones(7,1), 'k', ...
  'linewidth', 1.1)
hold on
theta = -2*pi:0.2:2*pi;     % cell edge
plot3(h*cos(theta), h*sin(theta), -0.5*ones(numel(theta),1), 'k--', ...
  'linewidth', 0.5)
theta = -2*pi:0.2:2*pi;     % cell nucleus edge
plot3(0.45*h*cos(theta), 0.45*h*sin(theta), -0.5*ones(numel(theta),1), ...
  'k--', 'linewidth', 0.5)
% ...visualize mesh above it
hold on
ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.65);
patch('Faces',R,'Vertices',V, 'FaceAlpha', 0.75, ...
  'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8]);
patch('Faces',R(ii,:), ...
'Vertices',V, 'FaceAlpha', 0.75,'FaceColor', ...
  graphics_color('bluish green'), 'EdgeColor',[0 0 0]);
dof = 28; %170;  % choose dof to emphasize and 'magnify'
vrtx = V(R(dof,:),:);
patch('Faces',R(dof,:), ...
'Vertices',V,'FaceColor',graphics_color('reddish purple'), ...
  'EdgeColor',[0 0 0]);
axis off
view([10 45 45])

% 'help' lines from smal voxel to big voxel
plot3([vrtx(3,1) +sqrt(3/4)*h], [vrtx(3,2) -h/2], [0 -0.5], ...
   'k--', 'linewidth', 1.1)
plot3([vrtx(6,1) -sqrt(3/4)*h], [vrtx(6,2) +h/2], [0 -0.5], ...
   'k--', 'linewidth', 1.1)
[X,Y] = meshgrid(-1:0.1:1);
% inner layer 'surface'
Z = -0.5*ones(size(X));
s = surf(X,Y,Z,'FaceAlpha',0.5, 'EdgeColor', 'none', 'FaceColor', ...
  [0.9 0.9 0.9]);

% edge and voxel center distance visualisation (e_ij, d_ij)
dof_2 = dof;     % choose two dofs on the side for clarity
vrtx_2 = V(R(dof_2,:),:);
dof_3 = dof-1;
vrtx_3 = V(R(dof_3,:),:);
% draw shared edge, e_ij
plot3([vrtx_2(6,1) vrtx_2(1,1)], [vrtx_2(6,2) vrtx_2(1,2)], [0 0], ...
  'r', 'linewidth', 2.5)
% draw distance between central nodes, d_ij
plot3([P(1,dof_2) P(1,dof_3)], [P(2,dof_2) P(2,dof_3)], [0 0], ...
  'b-|', 'linewidth', 2.5)

% clarify voxel v_i and v_j by text
text(P(1,dof_2)+0.18, P(2,dof_2)+0.16, 0.05, '$\Omega_i$', ...
  'Interpreter','latex')
text(P(1,dof_3)+0.03, P(2,dof_3)+0.16, 0.05, '$\Omega_j$', ...
  'Interpreter','latex')

% add text to figures
text(+0.8,-0.4,0.05, '$dt$', 'Interpreter','latex')
text(+0.8,-0.4,-0.45, '$d\tau$', 'Interpreter','latex')

% 'reactions here'...
text(+0.21,+0.15,-0.45, '$X \rightarrow Y$', 'Interpreter','latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 230 210]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

% tweak these parameters to fully remove figure whitespace...
t.Position = [-0.2 -0.05 1.4 1.09]; % left, bottom, right, top
