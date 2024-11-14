% HES1RED_GRID2D Reduced Hes1 ODE model over a fixed 2D grid.
%
%  See also HES1_GRID2D.
%
%  Reduction alternatives are as follows
%  
%  alternative 1: (x,y) = (M,n) [type 1]
%  alternative 2: (x,y) = (P,n) [type 1]
%  alternative 3: (x,y) = (M,N) [type 2]
%  alternative 4: (x,y) = (P,N) [type 2]
%  alternative 5: (x,y) = (M,D) [type 1]
%  alternative 6: (x,y) = (P,D) [type 1]
%  alternative 7: (x,y) = (M,P) [type 3]
%
%  Comment on notation in reduced model: 
%     (N,D) used here are equivalent to (x,y) in reference
%
%  Reference:
%     [1] G. Menz and S. Engblom: "Modelling Population-Level Hes1 Dynamics:
%         Insights from a Multi-Framework Approach", manuscript (2024).

% Gesina Menz 2024-02-23

if ~exist('save_figure', 'var')
  save_figure = 0;
end
if ~exist('save_data','var')
  save_data = 0;
end
if ~exist('alternative','var')
  alternative = 1;
end
if ~exist('par','var')
  par = hes1_params;
end
if ~exist('rand_seed', 'var')
  rand_seed = rng(1000)
else
  rng(rand_seed) % echo it
end

% set random seed
% rng(456)

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 20;
if mod(Nvoxels,2) == 1
  error("Nvoxels has to be an even number.")
end

% fetch discretization
[P,E,T,gradquotient] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% find neighbor matrix based on hexagonal mesh
N_op = dt_neighe(V,R);
N_op(N_op > 0) = 1; % set all non-zero entries to 1

% Hes1 DOFs over this mesh
Nmesh = size(P,2);
neigh = N_op*ones(Nmesh,1); % # of neighbors

% find layers on mesh
layers = mesh_layers(N_op,floor(Nvoxels/2+1));

% reduced Hes1 model on this grid
%   N' = F(N,D) - N
%   D' = v (G(N) - D)
%   where
%   F(x,y) = y / (a + x^k)
%   G(x) = 1 / (1 + b x^h)
%   N and D determined according to alternative chosen
%   e.g. for alternative 1: N <-> M, D <-> n

% reduction scalings for chosen alternative
% for alternative 1 (best fit): M0 is scale.x0, n0 is scale.y0
% for alternative 5: M0 is scale.x0, D0 is scale.y0
scale = hes1red_params(alternative,par);

F = @(x) 1 ./ (scale.a + x.^par.k);
G = @(x) 1./ (1 + scale.b * x.^par.h);
if alternative == 3 || alternative == 4
  % type 2 alternatives
  hes1 = @(N,D)[D .* F(N) - N; ...
      scale.v * (N_op*G(N)./neigh - D)];
elseif alternative == 7
  % type 3 alternative
  hes1 = @(N,D)[N_op*G(D)./neigh .* F(D) - N; ...
    scale.v * (N - D)];
else
  % type 1 alternatives
  hes1 = @(N,D)[((N_op*D)./neigh) .* F(N) - N; ...
    scale.v * (G(N) - D)];
end

% final ODE RHS
RHS = @(t,Y)hes1(Y(1:Nmesh),Y(Nmesh+1:end));

% solve
Tend = 6000/scale.tau;
tspan_ = linspace(0,Tend,401);
if ~exist('init','var') || numel(init) ~= 3*Nvoxels^2
  conc = hes1_conc;
  if alternative == 1 % x = M, y = n
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(3)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(5)/scale.y0;
  elseif alternative == 2 % x = P, y = n
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(4)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(5)/scale.y0;
  elseif alternative == 3 % x = M, y = N
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(3)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(2)/scale.y0;
  elseif alternative == 4 % x = P, y = N
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(4)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(2)/scale.y0;
  elseif alternative == 5 % x = M, y = D
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(3)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(1)/scale.y0;
  elseif alternative == 6 % x = P, y = D
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(3)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(1)/scale.y0;
  elseif alternative == 7 % x = M, y = P
    init = rand(2*Nmesh,1);
    init(1:Nmesh) = init(1:Nmesh).*conc(3)/scale.x0;
    init(Nmesh+1:2*Nmesh) = init(Nmesh+1:2*Nmesh).*conc(4)/scale.y0;
  else
    error('Choose an alternative between 1 and 7.')
  end
end
sol_ = ode15s(RHS,tspan_,init);
N_ = deval(sol_,tspan_,1:Nmesh);
D_ = deval(sol_,tspan_,Nmesh+1:2*Nmesh);

% sort data according to Hes1 protein at the end
[Nincr,idx] = sort(N_(:,end));
[~,ijmp] = max(diff(Nincr)); % cut at largest jump
idxlo = idx(1:ijmp); % index of low cells
idxhi = idx(ijmp+1:end); % index of high cells
idx_log = linspace(1,400,400);
idx_log(idxlo) = 0;
idx_log(idxhi) = 1;

% find 2 different groups for each consituent
D_lo = D_(idxlo,:);
D_hi = D_(idxhi,:);
N_lo = N_(idxlo,:);
N_hi = N_(idxhi,:);

if save_data == 1
  save('../data/hes1red_grid2D.mat','V','R','tspan_','N_',...
    'idx','idxlo','idxhi','scale','layers')
end

return;

% postprocessing
figure(1), clf,
plot(tspan_,mean(D_/scale.y0));
ylim([0 1]);
title('Mean Ngn2 per cell');
xlabel('time');

% appearance reduced system
figure(2),
max_ngn2 = max(D_,[],'all'); % find maximum Ngn2 level
division = [0,0.05,0.1,0.15,0.2,max_ngn2]/scale.y0;
cmap = colormap(parula(numel(division)-1));
axis_dim = [min(V),max(V)];
for i = 1:numel(tspan_)
  clf,
  patch('Faces',R,'Vertices',V, ...
    'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  patch('Faces',R(D_(:,i)/scale.y0 <= division(2),:), ...
    'Vertices',V,'FaceColor',cmap(2,:));
  for j = 3:numel(division)-1
    patch('Faces',R((D_(:,i)/scale.y0 > division(j-1) & D_(:,i)/scale.y0 <= division(j)),:), ...
      'Vertices',V, 'FaceColor',cmap(j-1,:));
  end
  patch('Faces',R(D_(:,i)/scale.y0 > division(numel(division)-1),:), ...
    'Vertices',V, 'FaceColor',cmap(numel(division)-1,:));
  axis([axis_dim(1) axis_dim(3) axis_dim(2) axis_dim(4)]);
  axis square, axis off
  title('Yellow: high Ngn2, Dark blue: low Ngn2');
  ticklabels=arrayfun(@(a)num2str(a),division,'uni',0);
  cb = colorbar('Ticks',linspace(0,1,numel(division)), ...
    'TickLabels',ticklabels);
  drawnow; pause(0.01);
end

% % the dynamics in a row
% figure(3), clf,
% ix = find(0 < P(2,:) & P(2,:) < 1e-1); % (line y ~ 0)
% %plot(P(1,ix),P(2,ix),'ko-');
% for i = 1:numel(tspan_)
%   clf,
%   stairs(P(1,ix),D_(ix,i),'k-'); hold on,
%   stairs(P(1,ix),N_(ix,i),'b--');
%   ylim([0 1]); yline(0.5,'r--');
%   drawnow; pause;
% end

% plot all Hes1 mRNA results to observe behaviour/final fate
figure(3), clf,
cmap = colormap(hsv(11));
for i = 1:11
  plot(tspan_,N_(layers{i,2},:)/scale.x0, 'Color',cmap(i,:))
  hold on
end
hold off
xlabel('time')
title('Hes1 mRNA behaviour in all cells')

% plot all Ngn2 results to observe behaviour/final fate
figure(4), clf,
cmap = colormap(hsv(11));
for i = 1:11
  plot(tspan_,D_(layers{i,2},:)/scale.y0, 'Color',cmap(i,:))
  hold on
end
hold off
xlabel('time')
title('Ngn2 behaviour in all cells')

% histogram at the end
figure(5), clf,
histogram(D_(:,end)/scale.y0);
title('Ngn2 behaviour at the end')

% histograms over time
figure(6), clf,
M = struct('cdata',{},'colormap',{});
for i=1:numel(tspan_)
  histogram(D_(:,i)/scale.y0,linspace(0,max_ngn2/scale.y0,100))
  title("Ngn2 levels at time " + tspan_(i))
  ylim([0 400])
  xticks(0:0.25:max_ngn2/scale.y0)
  drawnow; pause(0.01);
  M(i) = getframe(gcf);
end
if save_figure == 1
  movie2gif(M,{M([1:2 end]).cdata},'grid_ode_histograms.gif','delaytime',0.1,'loopcount',0);
end

% synchronisation for all cells
figure(7), clf,
plot(tspan_(1:150), abs(D_(:,1:150)/scale.y0-((N_op*D_(:,1:150)/scale.y0)./neigh)))
xlabel('time')
title('n-n_{mean} for all cells')

% mean synchronisation for all cells
figure(8), clf,
plot(tspan_, mean(abs(D_/scale.y0-((N_op*D_/scale.y0)./neigh))))
xlabel('time')
title('mean of n-n_{mean} over all cells')

% synchronisation over all cells
figure(9), clf,
plot(tspan_(1:150), abs(D_(:,1:150)/scale.y0-mean(D_(1:150)/scale.y0)))
xlabel('time')
title('n-mean(n) over all cells')