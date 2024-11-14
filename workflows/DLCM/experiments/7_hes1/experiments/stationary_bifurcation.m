%STATIONARY_BIFURCATION Hes1 ODE model over a 1D grid.
%   This experiment looks at the stationary solution(s) of the
%   Hes1-model at a 1D grid with periodic boundary conditions and as a
%   function of the parameters. By scaling the parameters, a
%   bifurcation diagram is obtained.

% S. Engblom 2024-09-24 (new model & parametrization)
% S. Engblom 2024-07-10

% cells live on a line of voxels
Nvoxels = 30;

% plots: hi+lo P(t)/plot stills/animate
output.hilo = 1;
output.stills = 0;
output.animate = 0;

% discretization: neighbor matrix
N_op = spdiags(ones(Nvoxels,1)*[1 0 1],-1:1,Nvoxels,Nvoxels);
N_op(1,[2 end]) = 1;
N_op(end,[1 end-1]) = 1;
N_op = N_op./sum(N_op,2); % average

% Hes1 model on this grid
ndof = 5;
par = hes1_params;
alpha = [par.alphaD par.alphaN par.alphaM par.alphaP par.alphan];
mu = [par.muD par.muN par.muM par.muP par.mun];
H = [par.KM par.Kn par.k par.h];

% determine critical bifurcation scaling
crit = @(s)l_criterion(alpha./[1 s 1 1 1],mu,H);
opts = optimset('Display','off');
[scrit,~,flag] = fzero(crit,500,opts);
if flag ~= 1, warning('Convergence issues.'); end
scrit

% Hes1 model with a parameter scaling of alphaN by factor s
[fun,funJ,funC] = hes1_buildODE;
RHS = @(s,t,Y)hes1_System(Y,alpha./[1 s 1 1 1],mu,H,fun,funC,N_op);

% "biased" initial data, then solve as an ODE
Y0 = repmat([1 1 1 1 1; 2 2 2 2 2]',1,Nvoxels/2)+1e-2*rand(ndof,Nvoxels);
tspan = linspace(0,5000,10);
[~,Y] = ode15s(@(t,y)RHS(1,t,y),tspan,Y0(:), ...
               odeset('RelTol',1e-8,'AbsTol',1e-8));

% continuation over parameter s
s = linspace(1,200);
s = [s(s < scrit) scrit s(s > scrit)];
Ystat = zeros(numel(Y0),numel(s));
Ystat(:,1) = Y(end,:)';
opts = optimset('Display','off','tolx',1e-8,'tolfun',1e-10);
for i = 2:numel(s)
  [Ystat(:,i),~,flag] = fsolve(@(y)RHS(s(i),0,y),Ystat(:,i-1),opts);
  if flag ~= 1, sprintf('Convergence issues, i = %d.',i); end
end

% solve by continuation backwards for the homogeneous stationary
Ystat_hom = zeros(numel(Y0),numel(s));
Ystat_hom(:,end) = Ystat(:,end);
for i = numel(s)-1:-1:1
  [Ystat_hom(:,i),~,flag] = ...
      fsolve(@(y)RHS(s(i),0,y),Ystat_hom(:,i+1),opts);
  if flag ~= 1, sprintf('Convergence issues, i = %d.',i); end
end

% investigate state P - normalized
Pstat = Ystat(4:ndof:end,:);
Pmax = max(Pstat(:));
Pstat_ = Pstat/Pmax;

Pstat_hom = Ystat_hom(4:ndof:end,:);
Pstat_hom_ = Pstat_hom/Pmax;

% bifurcation as a function of s
if output.hilo
  figure(1), clf,
  % find hi and lo from Pstat(:,1)
  [~,ix] = sort(Pstat(:,1));
  ixlo = ix(1:max(floor(0.25*Nvoxels),1));
  ixhi = ix(end-numel(ixlo)+1:end);
  jx = find(s <= scrit);
  kx = find(s >= scrit);
  semilogy(s(jx),mean(Pstat_(ixlo,jx),1),'b','LineWidth',2, ...
           'HandleVisibility','off');
  hold on,
  plot(s(jx),mean(Pstat_(ixhi,jx),1),'b','LineWidth',2);
  plot(s(jx),mean(Pstat_hom_(:,jx),1),'r--','LineWidth',2, ...
       'HandleVisibility','off');
  plot(s(kx),mean(Pstat_hom_(:,kx),1),'r','LineWidth',2);
  xline(scrit,'k--','HandleVisibility','off');

  % make this plot publishable
  set(gcf,'defaultTextInterpreter','latex');
  legend('Non-homogeneous','Homogeneous','Interpreter','latex', ...
         'location','SE');
  set(gca,'TickLabelInterpreter','latex');
  xlabel('$s$');
  ylabel('$P$');
  set(gca,'xtick',0:100:s(end),'ytick',0:1,'yticklabel',1,'fontsize',9);
  ax = axis; axis([s([1 end]) ax(3) 1]);
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[200 200 260 140]);
  %print -depsc ~/Desktop/fac_bifurcation.eps
end

% a few stills
if output.stills
  figure(2), clf,
  subplot(4,1,1);
  bar(1:Nvoxels,Pstat_(:,1));
  ylim([0 1]);
  subplot(4,1,2);
  bar(1:Nvoxels,Pstat_(:,20));
  ylim([0 1]);
  subplot(4,1,3);
  bar(1:Nvoxels,Pstat_(:,50));
  ylim([0 1]);
  subplot(4,1,4);
  bar(1:Nvoxels,Pstat_(:,100));
  ylim([0 1]);
  drawnow;
end

% animate
if output.animate
  figure(3), clf,
  for i = 1:numel(s)
    bar(1:Nvoxels,Pstat_(:,i));
    ylim([0 1]);
    drawnow;
    pause(0.1);
  end
end

% ----------------------------------------------------------------------
function res = l_criterion(alpha,mu,H)
%L_CRITERION Criterion for existence of non-homogeneous solutions etc.

[a,b] = hes1reduce(alpha,mu,H);
k = H(3); h = H(4);
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b.*x.^h);
df = @(x)-k*x.^(k-1).*f(x).^2;
dg = @(x)-h*b.*x.^(h-1).*g(x).^2;
x0 = hes1red_homogeneous(a,b,k,h);
res = f(x0).*dg(x0)-df(x0).*g(x0)+1;

end
% ----------------------------------------------------------------------
