%Morphogenesis 2D script file.
%   This file runs the Schnakenberg and the Brusselator models in
%   simple 2D geometries.
%
%   Reference:
%     [1] Y. Saygun: "Computational Stochastic Morphogenesis", MSc
%     thesis in Engineering Physics, Uppsala university (2015).

% S. Engblom 2017-03-08

%% (1) Schnakenberg

% build the geometry
C1 = [1 0 0 50]';
C2 = [1 0 0 15]';
gd = [C1 C2];
sf = 'C1-C2';
ns = char('C1','C2')';
G = decsg(gd,sf,ns);

% create the mesh
[P,E,T] = initmesh(G,'hmax',2.5);

% assemble the diffusion part
D_U = 1;
D_V = 40;
clear umod
umod = pde2urdme(P,T,{D_U D_V});
if exist('tspan','var')
  umod.tspan = tspan;
end
  
% not used
umod.sd = ceil(umod.sd);

if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf,
  pdegplot(G,'subdomainlabels','on'), axis equal

  figure(2), clf,
  pdemesh(P,E,T), axis tight, axis equal
end

umod = schnakenberg(umod);
umod.vol = 50/mean(umod.vol)*umod.vol;

% solve
umod = urdme(umod,'seed',17,'report',0);

% visualize using PDE Toolbox
umod = urdme2pde(umod);
if ~exist('plotting_off','var') || ~plotting_off
  figure(3), clf,
  pdesurf(umod.pde.P,umod.pde.T,umod.pde.U(1,:,end)');
  title('Schnakenberg: Concentration U');
  view(0,90), axis tight, axis square, colormap('parula')
  figure(4), clf,
  h = pdesurf(umod.pde.P,umod.pde.T,umod.pde.U(2,:,end)');
  title('Schnakenberg: Concentration V');
  view(0,90), axis tight, axis square, colormap('parula')
end

%% (2) Brusselator

% build the geometry
C1 = [1 0 0 25]';
gd = [C1];
sf = 'C1';
ns = char('C1')';
G = decsg(gd,sf,ns);

% create the mesh
[P,E,T] = initmesh(G,'hmax',2);

% assemble the diffusion part
D_U = 2;
D_V = 16;
clear vmod
vmod = pde2urdme(P,T,{D_U D_V});
if exist('tspan','var')
  vmod.tspan = tspan;
end

% not used
vmod.sd = ceil(vmod.sd);

vmod = brusselator(vmod);
vmod.vol = 100/mean(vmod.vol)*vmod.vol;

% solve
vmod = urdme(vmod,'seed',17,'report',0);

% visualize using PDE Toolbox
vmod = urdme2pde(vmod);
if ~exist('plotting_off','var') || ~plotting_off
  figure(5), clf,
  pdesurf(vmod.pde.P,vmod.pde.T,vmod.pde.U(1,:,end)');
  title('Brusselator: Concentration U');
  view(0,90), axis tight, axis square, colormap('parula')
  figure(6), clf,
  pdesurf(vmod.pde.P,vmod.pde.T,vmod.pde.U(2,:,end)');
  title('Brusselator: Concentration V');
  view(0,90), axis tight, axis square, colormap('parula')
end

return;

% possibly print to files

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
%print -depsc  ../../doc/fig/Schnakenberg_mesh.eps

figure(3),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
%print -depsc  ../../doc/fig/Schnakenberg.eps
