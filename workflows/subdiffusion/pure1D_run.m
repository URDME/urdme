% Subdiffusion 1D example.
%   Single species, pure subdiffusion in 1D.
%
%   Note: crowding/subdiffusion can be questioned as a concept in
%   1D. The purpose of this model is to test the methodology in a
%   simple case.

% S. Engblom 2017-03-24
% S. Engblom 2017-02-20 (Revision)
% S. Engblom 2014-05-09

clear umod;

% 1D interval [-1;1]
Nmesh = 50;
h = 2/Nmesh;

% diffusion matrix
e = ones(Nmesh,1);
D = spdiags([e -2*e e]/h^2,-1:1,Nmesh,Nmesh);
D(1,1) = -1/h^2;
D(end,end) = -1/h^2;
D = 0.01*D; % scaling of diffusion

% internal states approximation from pre-simulations
[gamma,A,p] = subdata(0.2);
Nstates = numel(gamma);

% initial- and equilibrium effective diffusion
[~,igammamax] = max(gamma);
gamma0 = gamma(igammamax);
gammabar = p'*gamma;

% given the internal state transfer matrix A, the diffusion
% coefficients gamma, the scalar diffusion operator D, compute the
% effective subdiffusion matrix
umod.D = subdmatrix(A,gamma,D);

% initial vector u0: all particles in the middle, and in the fastest
% diffusing state
umod.u0 = zeros(Nstates,Nmesh);
umod.u0(igammamax,Nmesh/2) = 1e2;

% unused
umod.N = sparse(Nstates,0);
umod.G = sparse(0,Nstates);
umod.vol = ones(1,Nmesh);
umod.sd = ones(1,Nmesh);

% output times
Tend = 20;
umod.tspan = linspace(0,Tend,101);

% solve
umod = urdme(umod);
vmod = umod;
vmod = urdme(vmod,'solver','uds');

% pure diffusion versions
wmod = umod;
wmod.u0 = sum(wmod.u0,1);
wmod.N = sparse(1,0);
wmod.G = sparse(0,1);
wmod.D = gamma0*D;
wmod1 = urdme(wmod,'solver','uds');
wmod.D = gammabar*D;
wmod2 = urdme(wmod,'solver','uds');

% plot regardless of inner state
msh = linspace(-1,1,Nmesh+1);
msh = 0.5*(msh(1:end-1)+msh(2:end));
t = umod.tspan;
U = reshape(umod.U,Nstates,Nmesh,[]);
V = reshape(vmod.U,Nstates,Nmesh,[]);
W1 = wmod1.U;
W2 = wmod2.U;

if exist('plotting_off','var') && plotting_off, return; end

figure(1), clf
M = struct('cdata',{},'colormap',{});
for i = 1:numel(t)
  hold off,
  plot(msh,W2(:,i)','k--'); hold on, ax = axis;
  plot(msh,W1(:,i)','k--');
  stairs(msh,sum(U(:,:,i),1)','b');
  plot(msh,sum(V(:,:,i),1)','r');
  axis([-0.5 0.5 0 ax(4)]), drawnow;
  M(i) = getframe;
end

% illustrative snapshots
j = 2;
for i = [6 51]
  figure(j), clf, j = j+1;
  plot(msh,W2(:,i)','k--'); hold on, ax = axis;
  plot(msh,W1(:,i)','k--');
  stairs(msh,sum(U(:,:,i),1)','b:');
  plot(msh,sum(V(:,:,i),1)','r');
  axis([-0.5 0.5 0 ax(4)]), drawnow;
end
figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 350 180 200]);
set(gca,'ytick',[]);

figure(3),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[300 350 180 200]);
set(gca,'ytick',[]);

% MSD-curves
figure(4), clf
msdV = (msh.^2)*squeeze(sum(V,1))./squeeze(sum(sum(V,1),2))';
msdW1 = (msh.^2)*W1./sum(W1,1);
msdW2 = (msh.^2)*W2./sum(W2,1);
loglog(vmod.tspan,msdV,'r','linewidth',2); hold on
loglog(wmod1.tspan,msdW1,'k--');
loglog(wmod2.tspan,msdW2,'k--');
exact = 2*0.01*t(end)*[gammabar gamma0];
axis([t(1) t(end) 0 exact(1)*2]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 400 150]);
set(gca,'xtick',[1e0 1e1],'ytick',[1e-2 1e-1],'yaxislocation','right');
h = xlabel('time');
set(h,'fontsize',8)
h = ylabel('MSD');
set(h,'fontsize',8)

% slope alpha
alpha = polyfit(log(t(2:end)),log(msdV(2:end)),1);
alpha = alpha(1)
tt = t(round([0.1 0.25]*numel(t)));
loglog(tt,exp(-6)*tt.^alpha,'b','linewidth',2);

return;

% print to file
figure(2),
print -depsc sample1d1.eps
figure(3),
print -depsc sample1d2.eps
figure(4),
print -depsc msd1d.eps
