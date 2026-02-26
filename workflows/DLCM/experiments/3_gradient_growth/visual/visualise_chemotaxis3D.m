%VISUALISE_CHEMOTAXIS Visualise chemotaxis model in 3D.

% E. Blom 2024-12-02

load chemotaxis_3D.mat

%% (1) Plot slices of chemotactic sensitivity
figure()
t = tiledlayout(3,2, 'Padding', 'none', 'TileSpacing', 'tight');

idx = find(strcmp([umod.solverargs{:}], 'mumod')) + 1; % find mumod
[X, Y, Z] = meshgrid(-1:0.01:1,-1:0.01:1, -1:0.01:1);

str = ['a)', 'b)', 'c)', 'd)'];
n = 1;  % loop counter

tscale = 3000; % 1/0.02*60 to get hours (1/0.02 rate to move 1 cell diam.)
tspan = [0 1]; % *10^5
for tt = [1 100]
  nexttile
  data = umod.solverargs{idx}{:}.U(1,:,tt);  % Chemical sensitivity
  UU = griddata(P(1,:), P(2,:), P(3,:), data, X, Y, Z);
  xslice = [-0.25,0,0.25]; yslice = []; zslice = [];
  s = slice(X,Y,Z,UU,xslice,yslice,zslice);
  s(1).EdgeColor = 'none';
  s(2).EdgeColor = 'none';
  s(3).EdgeColor = 'none';
  text(-0.3,+1,1.5, "$"+str(n:n+1)+"$", 'Interpreter','latex') % abcd)
  if tt == 1
    text(0.0,+1.2,1.5, "$t="+tspan(round(n/2))+"$ hours", ...
      'Interpreter','latex') % time
  else
    text(-0.0,+1.2,1.5, "$t="+round(tspan(round(n/2))*1e5/tscale)+"$ hours", ...
      'Interpreter','latex') % time
    colorbar('TickLabelInterpreter', 'Latex')
  end
  n = n+2;
  xlabel('x', 'Interpreter', 'Latex')
  ylabel('y', 'Interpreter', 'Latex')
  zlabel('z', 'Interpreter', 'Latex')
  set(gca,'TickLabelInterpreter',...
        'latex');
end

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
  colorbar off
  view([45,-45, 20])
  axis(0.8.*[-1.0 1.0, -1.0, 1.0, -1.0, 1.0])
  lightangle(0,30)
  text(-0.5,-1.25,1.2, "$"+str(n:n+1)+"$", 'Interpreter','latex')
  n = n+2;
  legend('', '', 'Single', 'Double', 'Interpreter', 'Latex')
end

nexttile([1 2])
% plot aggregation of sensitive cells in around origo
% get dofs centered around origo, enough to contain all sensitive cells
r0 = ((0.4^3*pi*4/3*1/4)/(4/3*pi))^(1/3);
r0_dof = (P(1,:).^2 + P(2,:).^2 + P(3,:).^2)<r0^2;
plot(umod.tspan/tscale, ...
  reshape(sum(sum(umod.solverargs{idx}{:}.U(:,r0_dof,:)>4.5))/176, ...
  [1 100]), 'LineWidth',2)
hold on
axis([0 max(umod.tspan/tscale), 0 1])
xlabel('time [hours]', 'Interpreter', 'Latex')
ylabel('$\bar{N}_s$', 'Interpreter', 'Latex')
yyaxis right
% plot nr of doubles
plot(umod.tspan/tscale,reshape(sum(umod.U(1,:,:)==2), ...
  [1 100]), 'linewidth', 2)
grid on
ylabel('$\sum_i(u_i = 2)$', 'Interpreter', 'Latex')
legend('Cell aggregation', 'Number of doubles', 'Interpreter', 'Latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 3*160]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
set(gca,'TickLabelInterpreter',...
        'latex');

% tweak these parameters to fully remove figure whitespace...
t.Position = [0.052+0.06 0.06 0.8983-0.13 0.9972-0.086]; % left, bottom, right, top

% uncomment to save:
%exportgraphics(t,'chemtax3D_cells.pdf')