% PANEL_PLOTS Panel plots and constituent plot included in Hes1 paper.

% G. Menz 2024-06-11

%% Configure fonts etc
fontsize = 9;
font = 'CMU Serif';

% times to show for both panel figures
% t reflects 60, 180, 270, 1800 minutes (1,2.75,4.5,30 hours)
t = [5,12,19,121];

% end time for panel figures
t_end = 141;

% time to consider before fate decision for calculating division 
% (about 900 minutes = 15 hours)
t_fd = 61;

% define constants for scaling of # molecules vs concentration
avogadro = 6.022*1e23; % for calculating mol
cell_vol = 50*1e-15; % cell volume of mouse embryonal stem cell (50 um^3)
conc = hes1_conc;
mean_molecules = 8104;

%% Panel figure Hes1@5 and Hes1@2

figure(1), clf,
tiled = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
figure_width = 441;
figure_height = 360;
set(gcf,'Position',[100 100 figure_width figure_height]);
set(gcf,'PaperPositionMode','auto','papersize',[figure_width ...
  figure_height]);

% row 1: plot grid behaviour of Hes1@5 at times t

% load Hes1@5 data if it exists, otherwise run simulation
if ~exist('../data/hes1_grid2D.mat','file')
  hes1_grid2D;
else
  load('../data/hes1_grid2D.mat')
end

% sort Hes1 mRNA into groups based on Hes1 protein
M_lo = M_(idxlo,:);
M_hi = M_(idxhi,:);

% save version with edges
M_edges = M_;
% get rid of edges and layer 2 (layer closest to edges)
M_([layers{1,2};layers{2,2};layers{3,2}],:) = [];
% indices for grid without corners, edges and layer 2
idx_layer1 = linspace(1,size(M_edges,1),size(M_edges,1));
idx_layer1 = idx_layer1(setdiff(idx_layer1,[layers{1,2};layers{2,2};...
  layers{3,2}]));
% save indices for original model without corners, edges and layer2
idxlo_orig = setdiff(idx_layer1,idxhi);
idxhi_orig = setdiff(idx_layer1,idxlo);

% rename tspan to not overwrite it later with reduced model
tspan_orig = tspan_;

% division line for high/low (average over first 900 minutes, i.e. 15 hours 
% so before fate decision)
division = mean(M_(:,1:t_fd),'all');

% t = 60
ax1 = nexttile;
set(ax1,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_edges(:,t(1)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_edges(:,t(1)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_(t(1))/60),'FontWeight','normal')

% t = 180
ax2 = nexttile;
set(ax2,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_edges(:,t(2)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_edges(:,t(2)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_(t(2))/60),'FontWeight','normal')

% t = 270
ax3 = nexttile;
set(ax3,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_edges(:,t(3)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_edges(:,t(3)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_(t(3))/60),'FontWeight','normal')

% t = 1800
ax4 = nexttile;
set(ax4,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_edges(:,t(4)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_edges(:,t(4)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_(t(4))/60),'FontWeight','normal')

% row 2: plot grid behaviour of Hes1@2 at times t

% load Hes1@2 data if it exists, otherwise run simulation
if ~exist('../data/hes1red_grid2D.mat','file')
  hes1red_grid2D;
else
  load('../data/hes1red_grid2D.mat')
end

% find 2 different groups for N
N_lo = N_(idxlo,:);
N_hi = N_(idxhi,:);

% save data with edges
N_edges = N_;
% get rid of edges and closest layer (layer 2)
N_([layers{1,2};layers{2,2};layers{3,2}],:) = [];
% find indices for model without corners, edges and layer2
idxlo = setdiff(idx_layer1,idxhi);
idxhi = setdiff(idx_layer1,idxlo);

% division line for high/low (average over first 300 minutes, i.e. before
% fate decision)
division = mean(N_(:,1:21),'all');

% t = 60
ax5 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(N_edges(:,t(1)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(N_edges(:,t(1)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 180
ax6 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(N_edges(:,t(2)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(N_edges(:,t(2)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 270
ax7 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(N_edges(:,t(3)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(N_edges(:,t(3)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 1800
ax8 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(N_edges(:,t(4)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(N_edges(:,t(4)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% row 3: plot behaviour over time for mean Hes1 protein

% find mean for high and low cells
M_top = M_edges(idxhi_orig,:);
M_bottom = M_edges(idxlo_orig,:);
N_top = N_edges(idxhi,:);
N_bottom = N_edges(idxlo,:);

mean_M_top = mean(M_top);
mean_M_bottom = mean(M_bottom);
mean_N_top = mean(N_top);
mean_N_bottom = mean(N_bottom);

ax9 = nexttile([1 4]);
plot(repmat(tspan_orig(t(1)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_orig(t(2)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_orig(t(3)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_orig(t(4)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
% plot Hes1@5
plot(tspan_orig,mean_M_top,'Color',graphics_color('blue'), ...
  'LineWidth',2)
hold on
plot(tspan_orig,mean_M_bottom,'Color',graphics_color('orange'), ...
  'LineWidth',2)
hold on
% plot Hes1@2
% remember: plot Hes1@2 with right scaling!
plot(tspan_orig,mean_N_top*scale.x0,'Color',graphics_color('blue'), ...
  'LineWidth',2,'LineStyle','--')
hold on
plot(tspan_orig,mean_N_bottom*scale.x0,'Color',graphics_color('orange'),...
  'LineWidth',2,'LineStyle','--')
set(ax9,'xlim',[0 tspan_orig(t_end)],'ylim',[0 0.1],'xtick',...
  0:180:tspan_orig(end),'xticklabel',(0:180:tspan_orig(end))/60,'Box',...
  'off','ytick',linspace(0,0.1,6),'FontName',font,'FontSize',fontsize-1)
xlabel('time [hrs]','FontName',font,'FontSize',fontsize)
ylabel(['Conc ', '[$\mu M$]'],'Interpreter','latex','FontName',font,...
  'FontSize',fontsize)

% to save
% exportgraphics(gcf,'../figures/ode_panel.pdf','Resolution',300)  

%% Panel figure URDME

% load or generate data for low volume
if ~exist('../data/hes1umod_vol1.mat','file')
  error('Run hes1umod2D_run.m with VOL = 1 and save the data.')
else
  % with random initial conditions
  load('../data/hes1umod_vol1_rand.mat')
end
Mspecies = size(umod.N,1);
VOL_low = VOL;
% calculate factor to multiply number of molecules with to get the soln in
% microM
factor_low = 1e6/avogadro./(VOL_low*1e-15);
M_low = umod.U(4:Mspecies:end,:).*factor_low;
P_low = umod.U(5:Mspecies:end,:).*factor_low;
% save version with edges
M_low_edges = M_low;
% disregard values at edges and the layer closest to the edges (layer 2)
M_low([layers{1,2};layers{2,2};layers{3,2}],:) = [];
tspan_low = umod.tspan;
% indices for grid without corners, edges and layer 2
idx_layer1 = linspace(1,size(M_low_edges,1),size(M_low_edges,1));
idx_layer1 = idx_layer1(setdiff(idx_layer1,[layers{1,2};layers{2,2};...
  layers{3,2}]));

% find high and low cells
% sort data according to Hes1 protein at the end
[Pincr,idx_low] = sort(P_low(:,end));
[~,ijmp] = max(diff(Pincr)); % cut at largest jump
idxlo_low = idx_low(1:ijmp); % index of low cells
idxhi_low = idx_low(ijmp+1:end); % index of high cells
% save indices for original model without corners, edges and layer2
idxlo_low = setdiff(idx_layer1,idxhi_low);
idxhi_low = setdiff(idx_layer1,idxlo_low);

% load data for high volume (cell size)
if ~exist('../data/hes1umod_vol50.mat','file')
  error('Run hes1umod2D_run.m with VOL = 50 and save the data.')
else
  % with random initial conditions
  load('../data/hes1umod_vol50_rand.mat')
end
VOL_high = VOL;
% calculate factor to multiply number of molecules with to get the soln in
% microM
factor_high = 1e6/avogadro./(VOL_high*1e-15);
M_high = umod.U(4:Mspecies:end,:).*factor_high;
P_high = umod.U(5:Mspecies:end,:).*factor_high;
% save version with edges
M_high_edges = M_high;
% disregard values at edges and the layer closest to the edges (layer 2)
M_high([layers{1,2};layers{2,2};layers{3,2}],:) = [];
tspan_high = umod.tspan;

% find high and low cells
% sort data according to Hes1 protein at the end
[Pincr,idx_high] = sort(P_high(:,end));
[~,ijmp] = max(diff(Pincr)); % cut at largest jump
idxlo_high = idx_high(1:ijmp); % index of low cells
idxhi_high = idx_high(ijmp+1:end); % index of high cells
% save indices for original model without corners, edges and layer2
idxlo_high = setdiff(idx_layer1,idxhi_high);
idxhi_high = setdiff(idx_layer1,idxlo_high);

figure(2), clf,
tiled = tiledlayout(3,4,'Padding','tight','TileSpacing','compact');
figure_width = 433;
figure_height = 360;
set(gcf,'Position',[100 100 figure_width figure_height]);
set(gcf,'PaperPositionMode','auto','papersize',[figure_width ...
  figure_height]);

% row 1: plot grid behaviour of URDME with low volume (VOL = 15) at times t

% division line for high/low (average over first 900 minutes, i.e. before
% fate decision)
division = mean(M_low(:,1:t_fd),'all');

% t = 60
ax1 = nexttile;
set(ax1,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_low_edges(:,t(1)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_low_edges(:,t(1)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(round(tspan_low(t(1))/60,2)),'FontWeight','normal')

% t = 180
ax2 = nexttile;
set(ax2,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_low_edges(:,t(2)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_low_edges(:,t(2)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_low(t(2))/60),'FontWeight','normal')

% t = 270
ax3 = nexttile;
set(ax3,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_low_edges(:,t(3)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_low_edges(:,t(3)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_low(t(3))/60),'FontWeight','normal')

% t = 1800
ax4 = nexttile;
set(ax4,'fontname',font,'fontsize',fontsize);
axis_dim = [min(V),max(V)];
patch('Faces',R(M_low_edges(:,t(4)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_low_edges(:,t(4)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off
title("t = "+num2str(tspan_low(t(4))/60),'FontWeight','normal')

% row 2: plot grid behaviour of URDME with high volume at times t

% division line for high/low (average over first 600 minutes, i.e. before
% fate decision)
division = mean(M_high(:,1:41),'all');

% t = 60
ax5 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(M_high_edges(:,t(1)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_high_edges(:,t(1)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 180
ax6 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(M_high_edges(:,t(2)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_high_edges(:,t(2)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 270
ax7 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(M_high_edges(:,t(3)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_high_edges(:,t(3)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% t = 1800
ax8 = nexttile;
axis_dim = [min(V),max(V)];
patch('Faces',R(M_high_edges(:,t(4)) <= division,:), ...
  'Vertices',V,'FaceColor',graphics_color('orange'));
hold on
patch('Faces',R(M_high_edges(:,t(4)) > division,:), ...
  'Vertices',V,'FaceColor',graphics_color('blue'));
axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
axis square
axis off

% row 3: plot behaviour over time for mean Hes1 protein

% find mean for high and low cells
Mlo_top = M_low_edges(idxhi_low,:);
Mlo_bottom = M_low_edges(idxlo_low,:);
Mhi_top = M_high_edges(idxhi_high,:);
Mhi_bottom = M_high_edges(idxlo_high,:);

mean_Mlo_top = mean(Mlo_top);
mean_Mlo_bottom = mean(Mlo_bottom);
mean_Mhi_top = mean(Mhi_top);
mean_Mhi_bottom = mean(Mhi_bottom);

ax9 = nexttile([1 4]);
plot(repmat(tspan_low(t(1)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_low(t(2)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_low(t(3)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(repmat(tspan_low(t(4)),2,1),[0 100],'Color',[0.8 0.8 0.8],...
  'LineStyle','--','LineWidth',1.5)
hold on
plot(tspan_low,mean_Mlo_top,'Color',...
  graphics_color('blue'),'LineStyle','--','LineWidth',2)
hold on
plot(tspan_low,mean_Mlo_bottom,'Color',...
  graphics_color('orange'),'LineStyle','--','LineWidth',2)
hold on
plot(tspan_low,mean_Mhi_top,'Color',...
  graphics_color('blue'),'LineWidth',2)
hold on
plot(tspan_low,mean_Mhi_bottom,'Color',...
  graphics_color('orange'),'LineWidth',2)
set(ax9,'xlim',[0 tspan_low(t_end)],'ylim',[0 0.12],'xtick',...
  0:180:tspan_low(end),'ytick',linspace(0,0.12,7),'xticklabel',...
  (0:180:tspan_low(end))/60,'Box','off','FontName',font,'Fontsize',...
  fontsize-1)
xlabel('time [hrs]','FontName',font,'Fontsize',fontsize)
ylabel('Conc $[\mu M]$','FontName',font,'Fontsize',fontsize,...
  'Interpreter','latex')

% to save
% exportgraphics(gcf,'../figures/urdme_panel.pdf','Resolution',300) 

%% Comparison of different constituents in Hes1@5

figure(3), clf,
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 350 180]);

% load Hes1@5 data if it exists, otherwise run simulation
if ~exist('../data/hes1_grid2D.mat','file')
  hes1_grid2D;
else
  load('../data/hes1_grid2D.mat')
end
meanD = mean(D_);
meanN = mean(N_);
meanM = mean(M_);
meanP = mean(P_);
meann = mean(n_);

% plot mean of all constituents up to time 1000 minutes (16.67 hrs)
cut1 = tspan_(12); cut2 = tspan_(14); cut3 = tspan_(27);
ycut = [1e-10 max(mean(N_))];
semilogy([cut1 cut1],ycut,'Color',[0.8 0.8 0.8],'LineStyle','--',...
  'LineWidth',1.5,'HandleVisibility','off');
hold on
semilogy([cut2 cut2],ycut,'Color',[0.8 0.8 0.8],'LineStyle','--',...
  'LineWidth',1.5,'HandleVisibility','off');
hold on
semilogy([cut3 cut3],ycut,'Color',[0.8 0.8 0.8],'LineStyle',':',...
  'LineWidth',1.5,'HandleVisibility','off');
hold on
semilogy(tspan_,meanD,'Color',graphics_color('bluish green'),'Linewidth',...
  2,'DisplayName','Dll1')
hold on
semilogy(tspan_,meanN,'Color',graphics_color('orange'),'Linewidth',...
  2,'DisplayName','Notch')
hold on
semilogy(tspan_,meanM,'Color',graphics_color('reddish purple'),'Linewidth',...
  2,'DisplayName','Hes1 mRNA')
hold on
semilogy(tspan_,meanP,'Color',graphics_color('sky blue'),'Linewidth',...
  2,'DisplayName','Hes1 protein')
hold on
semilogy(tspan_,meann,'Color',graphics_color('vermillion'),'Linewidth',...
  2,'DisplayName','Ngn2')
hold on
semilogy(cut1,meanM(12),'ro','MarkerSize',8,'Linewidth',1.5,...
  'HandleVisibility','off')
hold on
semilogy(cut2,meanP(14),'ro','MarkerSize',8,'Linewidth',1.5,...
  'HandleVisibility','off')
hold on
semilogy(cut3,meanP(27),'kx','MarkerSize',8,'Linewidth',2,...
  'HandleVisibility','off')
hold on
semilogy(cut3,meanD(27),'kx','MarkerSize',8,'Linewidth',2,...
  'HandleVisibility','off')
hold on
semilogy(cut3,meann(27),'kx','MarkerSize',8,'Linewidth',2,...
  'HandleVisibility','off')
legend('Location','eastoutside','FontName',font,'fontsize',fontsize-1)
set(gca,'xtick',0:180:tspan_(68),'xticklabel',(0:180:tspan_(68))/60,...
  'FontName',font,'fontsize',fontsize-1)
xlim([0 1000])
ylim([0.000302341 0.8])
xlabel('time [hrs]','FontName',font,'FontSize',fontsize)
ylabel('Conc $[\mu M]$','FontName',font,'FontSize',fontsize,...
  'Interpreter','latex')

% to save
% exportgraphics(gca,'../figures/ode_constituents_comparison.pdf','Resolution',300) 