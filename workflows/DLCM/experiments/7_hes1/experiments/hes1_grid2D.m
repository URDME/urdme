% HES1_GRID2D Full Hes1 ODE model over a fixed 2D grid.

% Gesina Menz 2023-10-02

if ~exist('save_figure', 'var')
  save_figure = 0;
end
if ~exist('save_data','var')
  save_data = 0;
end
if ~exist('par','var')
  par = hes1_params;
end
if ~exist('rand_seed', 'var')
  rand_seed = rng(1000);
else
  rng(rand_seed) % echo it
end

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 20;
if mod(Nvoxels,2) == 1
  error("Nvoxels has to be an even number.")
end

% fetch discretization
[P,E,T,gradquotient] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% find coupling and neighbor matrix based on hexagonal mesh
W = dt_neighe(V,R);
N_op = W;
N_op(N_op > 0) = 1; % set all non-zero entries to 1

% Hes1 DOFs over this mesh
Nmesh = size(P,2);
neigh = N_op*ones(Nmesh,1); % # of neighbors

% find layers on mesh
layers = mesh_layers(N_op,floor(Nvoxels/2+1));

% Hes1 model on this grid
[fun,funJ,funC] = hes1_buildODE; % prepare once

% solve
if ~exist('init','var') || numel(init) ~= 5*Nvoxels^2
  % random initial values scaled to wanted concentrations to achieve
  % proper behaviour
  conc = hes1_conc;
  init = reshape(rand(5,Nmesh).*conc,[],1);
end
Tend = 6000;
tspan_ = linspace(0,Tend,401);
alpha = [par.alphaD par.alphaN par.alphaM par.alphaP par.alphan];
mu = [par.muD par.muN par.muM par.muP par.mun];
H = [par.KM par.Kn par.k par.h];

% optional use of Jacobian:
opts = odeset('AbsTol',1e-6,'RelTol',1e-8, ...
  'Jacobian',@(t,y)hes1_Jacobian(y,alpha,mu,H,funJ,funC,W));

% final solve
[~,sol_] = ode15s(@(t,y)hes1_System(y,alpha,mu,H, ...
  fun,funC,W),tspan_,init(:),opts);
D_ = sol_(:,1:5:end).';
N_ = sol_(:,2:5:end).';
M_ = sol_(:,3:5:end).';
P_ = sol_(:,4:5:end).';
n_ = sol_(:,5:5:end).';

% sort data according to Hes1 protein at the end
[Pincr,idx] = sort(P_(:,end));
[~,ijmp] = max(diff(Pincr)); % cut at largest jump
idxlo = idx(1:ijmp); % index of low cells
idxhi = idx(ijmp+1:end); % index of high cells
idx_log = linspace(1,400,400);
idx_log(idxlo) = 0;
idx_log(idxhi) = 1;

% sort each constituent into groups based on Hes1 protein
D_lo = D_(idxlo,:);
D_hi = D_(idxhi,:);
N_lo = N_(idxlo,:);
N_hi = N_(idxhi,:);
M_lo = M_(idxlo,:);
M_hi = M_(idxhi,:);
P_lo = P_(idxlo,:);
P_hi = P_(idxhi,:);
n_lo = n_(idxlo,:);
n_hi = n_(idxhi,:);

% find entropy over time
prob_dist = zeros(size(n_,1),5);
prob_dist_clust1 = zeros(size(N_lo,1),5);
prob_dist_clust2 = zeros(size(N_hi,1),5);
entropy = zeros(numel(tspan_),15);
for i = 1:numel(tspan_)
  for j = 1:size(n_,1) % over all cells
    prob_dist(j,1) = D_(j,i)/sum(D_(:,i));
    prob_dist(j,2) = N_(j,i)/sum(N_(:,i));
    prob_dist(j,3) = M_(j,i)/sum(M_(:,i));
    prob_dist(j,4) = P_(j,i)/sum(P_(:,i));
    prob_dist(j,5) = n_(j,i)/sum(n_(:,i));
  end
  for j = 1:size(N_lo,1) % over low Ngn2 cells
    prob_dist_clust1(j,1) = D_lo(j,i)/sum(D_lo(:,i));
    prob_dist_clust1(j,2) = N_lo(j,i)/sum(N_lo(:,i));
    prob_dist_clust1(j,3) = M_lo(j,i)/sum(M_lo(:,i));
    prob_dist_clust1(j,4) = P_lo(j,i)/sum(P_lo(:,i));
    prob_dist_clust1(j,5) = n_lo(j,i)/sum(n_lo(:,i));
  end
  for j = 1:size(N_hi,1) % over high Ngn2 cells
    prob_dist_clust2(j,1) = D_hi(j,i)/sum(D_hi(:,i));
    prob_dist_clust2(j,2) = N_hi(j,i)/sum(N_hi(:,i));
    prob_dist_clust2(j,3) = M_hi(j,i)/sum(M_hi(:,i));
    prob_dist_clust2(j,4) = P_hi(j,i)/sum(P_hi(:,i));
    prob_dist_clust2(j,5) = n_hi(j,i)/sum(n_hi(:,i));
  end
  entropy(i,1:5) = - sum(prob_dist.*log(prob_dist));
  entropy(i,6:10) = - sum(prob_dist_clust1.*log(prob_dist_clust1));
  entropy(i,11:15) = - sum(prob_dist_clust2.*log(prob_dist_clust2));
end

if save_data == 1
  save('../data/hes1_grid2D.mat','V','R','tspan_','D_','N_','M_','P_',...
    'n_','idx','idxlo','idxhi','layers','N_op')
end

% return;

% find power spectrum for two groups
% power_spectrum1 = pspectrum(P_(idx == 1));
% power_spectrum2 = pspectrum(P_(idx == 2));
% figure(1), clf,
% pspectrum(P_(idx == 1))
% figure(2), clf,
% pspectrum(P_(idx == 2));
% 
% % postprocessing
% figure(1), clf,
% plot(tspan_,mean(n_));
% ylim([0 1]);
% title('Mean Ngn2 per cell');
% xlabel('time');

% find mean expression for all constituents
meanD = mean(D_);
meanN = mean(N_);
meanM = mean(M_);
meanP = mean(P_);
meann = mean(n_);

% calculate period of mean oscillations by finding difference between local
% maxima (between 2nd and 3rd maximum to avoid issues with first oscillation
% being different due to initial condition)
% should be between 120-180 as periods of oscillations for Hes1 protein
[y1,t1] = findpeaks(meanP(1:end/2),tspan_(1:end/2));
fprintf('Measured period: %2.1fh \n',(t1(3)-t1(2))/60);

figure(1),clf,
% plot mean of all constituents up to time 1000 minutes (16.67 hrs)
% plot(repmat(tspan_(24),2,1),[0 100],'Color',[0.8 0.8 0.8],'LineStyle','--',...
%   'LineWidth',1.5,'HandleVisibility','off');
% hold on
% plot(repmat(tspan_(26),2,1),[0 100],'Color',[0.8 0.8 0.8],'LineStyle','--',...
%   'LineWidth',1.5,'HandleVisibility','off');
% hold on
semilogy(tspan_,mean(D_),'Color',graphics_color('bluish green'),'Linewidth',...
  2,'DisplayName','Dll1')
hold on
semilogy(tspan_,mean(N_),'Color',graphics_color('orange'),'Linewidth',...
  2,'DisplayName','Notch')
hold on
semilogy(tspan_,meanM,'Color',graphics_color('reddish purple'),'Linewidth',...
  2,'DisplayName','Hes1 mRNA')
hold on
semilogy(tspan_,meanP,'Color',graphics_color('sky blue'),'Linewidth',...
  2,'DisplayName','Hes1 protein')
hold on
semilogy(tspan_,mean(n_),'Color',graphics_color('vermillion'),'Linewidth',...
  2,'DisplayName','Ngn2')
hold on
plot(t1(2:3),repmat(y1(2),2,1),'b--.', 'LineWidth',2,'MarkerSize',10,...
  'DisplayName','measured period');
% hold on
% plot(tspan_(24),meanM(24),'ro','MarkerSize',8,'Linewidth',1.5,...
%   'HandleVisibility','off')
% hold on
% plot(tspan_(26),meanP(26),'ro','MarkerSize',8,'Linewidth',1.5,...
%   'HandleVisibility','off')
legend('Location','southeast')
set(gca,'xtick',0:240:tspan_(end),'xticklabel',(0:240:tspan_(end))/60)
% ylim([0 5])
% xlim([0 2500])
xlabel('time [hrs]')
ylabel('Conc [$\mu M$]','Interpreter','latex')

return;

% appearance
figure(2),clf,
max_ngn2 = max(n_,[],'all'); % find maximum Ngn2 level
% division = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.45,max_ngn2];
division = mean(M_(:,1:141),'all');
% cmap = colormap(parula(numel(division)-1));
axis_dim = [min(V),max(V)];
for i = 1:numel(tspan_)
  clf,
  patch('Faces',R,'Vertices',V, ...
    'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
%   patch('Faces',R(n_(:,i) <= division(2),:), ...
%     'Vertices',V,'FaceColor',cmap(2,:));
%   for j = 3:numel(division)-1
%     patch('Faces',R((n_(:,i) > division(j-1) & n_(:,i) <= division(j)),:), ...
%       'Vertices',V, 'FaceColor',cmap(j-1,:));
%   end
%   patch('Faces',R(n_(:,i) > division(numel(division)-1),:), ...
%     'Vertices',V, 'FaceColor',cmap(numel(division)-1,:));
  patch('Faces',R(n_(:,i) <= division,:), ...
    'Vertices',V,'FaceColor',graphics_color('blue'));
  hold on
  patch('Faces',R(n_(:,i) > division,:), ...
    'Vertices',V, 'FaceColor',graphics_color('yellow'));
  axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
  axis square, axis off
  title('Yellow: high Ngn2, Dark blue: low Ngn2');
  ticklabels=arrayfun(@(a)num2str(a),division,'uni',0);
  cb = colorbar('Ticks',linspace(0,1,numel(division)), ...
    'TickLabels',ticklabels);
  drawnow; pause(0.01);
end

% entropy over time
figure(3), clf,
subplot(3,1,1); % in all cells
plot(tspan_, entropy(:,1:2))
hold on
plot(tspan_, entropy(:,3))
hold on
plot(tspan_, entropy(:,4), '--')
hold on
plot(tspan_, entropy(:,5), '--')
title('Entropy over all cells')
legend({'Dll1', 'Notch', 'Hes1 mRNA', 'Hes1 protein', 'Ngn2'})

subplot(3,1,2); % in low Ngn2 cells
plot(tspan_, entropy(:,6:7))
hold on
plot(tspan_, entropy(:,8))
hold on
plot(tspan_, entropy(:,9), '--')
hold on
plot(tspan_, entropy(:,10), '--')
title('Entropy over low Ngn2 cells')

subplot(3,1,3); % in high Ngn2 cells
plot(tspan_, entropy(:,11:12))
hold on
plot(tspan_, entropy(:,13))
hold on
plot(tspan_, entropy(:,14), '--')
hold on
plot(tspan_, entropy(:,15), '--')
title('Entropy over high Ngn2 cells')

% plot all Hes1 results to observe behaviour/final fate
figure(4), clf,
cmap = colormap(hsv((Nvoxels/2)+1));
for i = 1:(Nvoxels/2)+1
  plot(tspan_,P_(layers{i,2},:), 'Color',cmap(i,:))
  hold on
end
hold off
xlabel('time [h]')
xticks(0:300:tspan_(end))
xticklabels((0:300:tspan_(end))/60)
title('Hes1 behaviour in all cells')

% histogram at the end
figure(5), clf,
histogram(n_(:,end));
title('Ngn2 behaviour at the end')

% histograms over time
figure(6), clf,
M = struct('cdata',{},'colormap',{});
for i=1:numel(tspan_)
  histogram(n_(:,i),linspace(0,max_ngn2,100))
  title("Ngn2 levels at time " + tspan_(i))
  ylim([0 400])
  xticks(0:0.25:max_ngn2)
  drawnow; pause(0.01);
  M(i) = getframe(gcf);
end
if save_figure == 1
  movie2gif(M,{M([1:2 end]).cdata},'grid_ode_histograms.gif','delaytime',0.1,'loopcount',0);
end

% synchronisation for all cells
figure(7), clf,
plot(tspan_(1:150), abs(D_(:,1:150)-((N_op*D_(:,1:150))./neigh)))
xlabel('time')
title('D-D_{mean} for all cells')

% mean synchronisation for all cells
figure(8), clf,
plot(tspan_, mean(abs(D_-((N_op*D_)./neigh))))
xlabel('time')
title('mean of D-D_{mean} over all cells')

% synchronisation over all cells
figure(9), clf,
plot(tspan_(1:150), abs(D_(:,1:150)-mean(D_(1:150))))
xlabel('time')
title('D-mean(D) over all cells')

% plot mean Hes1 results to observe behaviour/final fate
figure(10), clf,
plot(tspan_,mean(P_(idx == 1,:)), 'r')
hold on
plot(tspan_,mean(P_(idx == 2,:)), 'b')
hold off
xlabel('time')
title('Mean Hes1 behaviour over all cells')