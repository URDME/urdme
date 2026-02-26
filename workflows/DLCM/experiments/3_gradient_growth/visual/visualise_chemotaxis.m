%VISUALISE_CHEMOTAXIS Visualise samples of various chemotaxis models.

% E. Blom 2024-12-03

%% (1) Plot slices of chemotactic sensitivity
figure()
t = tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'none');
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
cb = colorbar('location', 'westoutside', 'TickLabelInterpreter',...
        'latex');
set(cb,'Position',[0.05 0.65 .015 0.33])

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

nexttile([1 4]) % plot quantitative results below

% visualise movement in x of center of mass
for n = 1:4
  load("chem" +n+ "_mean.mat")
  % $$$ use to extract position mean and std from a simulation
  % $$$ for tt = 1:numel(umod.tspan)
  % $$$   adof = find(umod.U(1,:,tt) > 0);
  % $$$   xymean(tt,:) = mean(umod.U(1,adof,tt).*[P(1,adof); P(2,adof)],2);
  % $$$ end

  % unit scaling
  % estimate of speed due to chemotaxis only: 0.05*0.0115*0.02
  % where 0.05*0.0115 is speed to move one cell diameter (grad(Q)*e_ij).
  % we do rough units by setting this rate to be one minute.
  % => 1740 time steps = 1 minute.
  tscale = 1740*60; % one hour per tscale time steps
  xscale = 2;       % x = 1 => 1/2 millimeter

  % plot
  %yyaxis left
  if n == 4
    tscale = tscale/0.0150; % scale 4th model to scale of 1-3rd
  end
  plot(tspan/tscale, -xymean(:,1)/xscale, 'LineWidth', 2)
  hold on;

  if n == 4
    tscale = tscale*0.0150; % use OG scale for reference line
  end
end

%yyaxis left
plot(tspan/tscale, 1.15e-5*tspan/xscale, 'k--', 'LineWidth', 1.5)
xlabel('time [hours]', 'Interpreter', 'Latex')
ylabel('distance in x [mm]', 'Interpreter', 'Latex')
grid on

legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Reference', ...
  'Interpreter', 'Latex')
axis([0 0.28 0 0.23])

% post-process figure
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 400 250]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
set(gca,'TickLabelInterpreter',...
        'latex');

% tweak these parameters to fully remove figure whitespace...
t.Position = [0.095 0.12 0.9 0.93]; % left, bottom, right, top
%t.Position = [0 -0.15 0.885 1];

% uncomment to save:
%exportgraphics(t,'chemtax_variety.pdf')
