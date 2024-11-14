% HES1UMOD_BIFURCATION Find bifurcation plot for RDME model.

% S. Engblom 2024-10-26 (Revision)
% G. Menz 2024-10-18

if ~exist('save_data','var')
  save_data = false;
end

% scaling s of alphaN; use svec from stationary_3cell
load ../data/stationary_3cell.mat
s = [svec svec(end:-1:1)]; % "fwd-bwd"

% max #iterations and tolerance for stationary
MAXITER = 50;
TOL = 0.01;

% build the model using external script
solve = 0;
VOL = 20;
% 3-cell problem, small, or large case
%Nvoxels = '3';
%Nvoxels = 5;
Nvoxels = 20;
hes1umod2D_run;

% modify the struct a bit:
umod.solve = 1; vmod.solve = umod.solve;
umod.parse = 1; vmod.parse = umod.parse;
umod.compile = 0; vmod.compile = umod.compile;
umod.report = 0; vmod.report = umod.report;
umod.seed = 241018;
tspan = linspace(0,24*60,10); % initial larger tspan
umod.tspan = tspan; vmod.tspan = umod.tspan;

% D will be modified using the (unscaled) D0
D0 = umod.D*umod.gdata;

% output approximately stationary P
saveP = zeros(size(umod.u0,2),numel(s));
savePuds = zeros(size(umod.u0,2),numel(s));

% solve model sequentially with scaling s
for i = 1:numel(s)
  disp(['Solving with scaling = ' num2str(s(i)) '... ']);

  % inject new scaling
  umod.seed = umod.seed+1;
  umod.gdata = s(i); vmod.gdata = umod.gdata;
  umod.D = D0/s(i); vmod.D = umod.D;

  % "shake" the initial data a bit
  umod.u0 = umod.u0.*exp(TOL*randn(size(umod.u0)));
  vmod.u0 = vmod.u0.*exp(TOL*randn(size(vmod.u0)));

  % solve until stationary
  for k = 1:MAXITER
    umod.seed = umod.seed+1;
    umod = urdme(umod);
    vmod = urdme(vmod);
    % note: criterion placed on vmod/uds, not umod/nsm
    if norm(vmod.u0(:)-vmod.U(:,end)) < TOL*norm(vmod.U(:,end))
      break;
    end
    umod.u0 = reshape(umod.U(:,end),size(umod.u0));
    vmod.u0 = reshape(vmod.U(:,end),size(vmod.u0));
  end
  if k == MAXITER
    warning('Stationary tolerance not fulfilled');
  end

  % save P
  saveP(:,i) = umod.U(5:Mspecies:end,end);
  savePuds(:,i) = vmod.U(5:Mspecies:end,end);
  
  % prepare for next round
  umod.u0 = reshape(umod.U(:,end),size(umod.u0)); 
  vmod.u0 = reshape(vmod.U(:,end),size(vmod.u0));
  umod.tspan = [0 4*60]; vmod.tspan = umod.tspan;

  disp('   ...done.');
end

% postprocess
l_diam = @(X)max(X,[],1)-min(X,[],1);

% since it is less noisy, base the sorting on the UDS-result:
P_ = sort(savePuds);
ixlo = 1:ceil(0.25*size(P_,1));
ixhi = size(P_,1)-numel(ixlo)+1:size(P_,1);
meanPuds = [mean(P_(ixlo,:),1); mean(P_(ixhi,:),1)];
diamPuds = [l_diam(P_(ixlo,:)); l_diam(P_(ixhi,:))];
rngPuds = [P_(ixlo([1 end]),:); P_(ixhi([1 end]),:)];
% judge overlap of P(ixlo) and P(ixhi):
ix = meanPuds(1,:)+1.5*diamPuds(1,:) >= meanPuds(2,:)-1.5*diamPuds(2,:);
% (another variant:)
%ix = rngPuds(2,:)+diamPuds(1,:) >= rngPuds(3,:)-diamPuds(2,:);
% logical negation of index used for plotting, but non-empty nix
% overlaps ix by 1 for visual purposes:
nix = ~ix;
if any(nix)
  nix(find(ix,1,'first')) = 1;
  nix(find(ix,1,'last')) = 1;
end
% transform into indices:
nix = find(nix);
ix = find(ix);
% overlapping part where there seems to be only one solution:
MeanPuds = nan(1,size(meanPuds,2));
DiamPuds = nan(1,size(meanPuds,2));
RngPuds = nan(2,size(meanPuds,2));
MeanPuds(:,ix) = mean(P_(:,ix));
DiamPuds(:,ix) = l_diam(P_(:,ix));
RngPuds(:,ix) = [min(P_(:,ix),[],1); max(P_(:,ix),[],1)];

P_ = sort(saveP);
meanP = [mean(P_(ixlo,:),1); mean(P_(ixhi,:),1)];
diamP = [l_diam(P_(ixlo,:)); l_diam(P_(ixhi,:))];
rngP = [P_(ixlo([1 end]),:); P_(ixhi([1 end]),:)];
MeanP = nan(1,size(meanP,2));
DiamP = nan(1,size(meanP,2));
RngP = nan(2,size(meanP,2));
MeanP(:,ix) = mean(P_(:,ix));
DiamP(:,ix) = l_diam(P_(:,ix));
RngP(:,ix) = [min(P_(:,ix),[],1); max(P_(:,ix),[],1)];

if save_data
  disp('Saving...')
  save("../data/hes1_bifurcation_Nvoxels_"+ ...
       num2str(Nvoxels)+"_VOL_"+num2str(VOL)+".mat", ...
       's','saveP','savePuds', ...
       'meanP','diamP','rngP','meanPuds','diamPuds','rngPuds', ...
       'MeanP','DiamP','RngP','MeanPuds','DiamPuds','RngPuds', ...
       'ix','nix', ...
       'scrit','VOL');
  disp('   ...saved.')
end

% visualization i = "fwd-bwd"
for i = 1:2:3
  figure(i), clf,
  jx = nix;
  kx = ix;
  if i == 1
    jx = jx(jx <= numel(s)/2);
    kx = kx(kx <= numel(s)/2);
  else
    jx = jx(jx > numel(s)/2);    
    kx = kx(kx > numel(s)/2);
  end
  % NSM
  P0 = max(meanP(:));
  h = plot(s(jx),meanP(1,jx)/P0,'r', ...
           'LineWidth',2,'HandleVisibility','off');hold on,
  if ~isempty(h)
    col = get(h,'Color');
    errorshade(s(jx),max(rngP(1,jx)/P0,1e-4),max(rngP(2,jx)/P0,1e-4),col);
    errorshade(s(jx),max(rngP(3,jx)/P0,1e-4),max(rngP(4,jx)/P0,1e-4),col);
  end
  % plot again or else the pdf does not work (transparency):
  plot(s(jx),meanP(1,jx)/P0,'r', ...
           'LineWidth',2);
  plot(s(jx),meanP(2,jx)/P0,'r','LineWidth',2,'HandleVisibility','off');
  h = plot(s(kx),MeanP(kx)/P0,'r','LineWidth',2,'HandleVisibility','off');
  hold on,
  if ~isempty(h)
    col = get(h,'Color');
    errorshade(s(kx),RngP(1,kx)/P0,RngP(2,kx)/P0,col);
  end
  plot(s(kx),MeanP(kx)/P0,'r','LineWidth',2,'HandleVisibility','off');
  set(gcf,'defaultTextInterpreter','latex');
  xlabel('$s$');
  ylabel('$P$');
  xline(scrit,'k','HandleVisibility','off');
  xline(scrit2,'k-.','HandleVisibility','off');
  axis([min(s) max(s) 1e-3 1]);
  set(gca,'defaultTextInterpreter','latex');
  set(gca,'xtick',0:50:max(s),'ytick',0:1,'yticklabel',0:1,'fontsize',9);

  % UDS
  P0 = max(meanPuds(:));
  figure(i+1), clf
  h = plot(s(jx),meanPuds(1,jx)/P0,'b','LineWidth',2); hold on,
  plot(s(jx),meanPuds(2,jx)/P0,'b','LineWidth',2,'HandleVisibility','off');
  if ~isempty(h)
    col = get(h,'Color');
    errorshade(s(jx),max(rngPuds(1,jx)/P0,1e-4),max(rngPuds(2,jx)/P0,1e-4),col);
    errorshade(s(jx),max(rngPuds(3,jx)/P0,1e-4),max(rngPuds(4,jx)/P0,1e-4),col);
  end
  h = plot(s(kx),MeanPuds(kx)/P0,'b','LineWidth',2,'HandleVisibility','off');
  if ~isempty(h)
    errorshade(s(kx),RngPuds(1,kx)/P0,RngPuds(2,kx)/P0,col);
  end
  set(gcf,'defaultTextInterpreter','latex');
  xlabel('$s$');
  ylabel('$P$');
  xline(scrit,'k','HandleVisibility','off');
  xline(scrit2,'k-.','HandleVisibility','off');
  axis([min(s) max(s) 0 1]);
  set(gca,'defaultTextInterpreter','latex');
  set(gca,'xtick',0:50:max(s),'fontsize',9);
end

% continue with Fig. 3 a bit:
figure(3),
set(gca,'YScale','log'); % for consistenty with other plots
axis([min(s) max(s) 1e-3 1]);
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',0:50:max(s), ...
        'ytick',[1e-3 1e-2 1e-1 1], ...
        'yticklabel',{'10^{-3}' '10^{-2}' '10^{-1}' '1'},'fontsize',9);
ix_ = find(lam(2,:) <= 0);
plot(svec(ix_),X1(1,ix_),'b','LineWidth',1);
plot(svec(ix_),X1(2,ix_),'b', ...
     'LineWidth',1);
ix_ = find(lam(1,:) <= 0);
plot(svec(ix_),X0(ix_),'b','LineWidth',1,'HandleVisibility','off');

% make plot publishable:
legend('RDME','ODE/$W_3$','location','SE','Interpreter','latex');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[200 200 260 160]);
%print -depsc ~/Desktop/rdme_bifurcation.eps
