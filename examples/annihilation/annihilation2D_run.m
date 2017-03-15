% Annihilation 2D example.
%   Species A and B are created in different subdomains and annihilate
%   when they meet. This example shows how one can set up a simple 2D
%   simulation using the URDME interface to PDE Toolbox.

% S. Engblom 2017-02-21

% generate URDME model
clear umod
[umod,G,P,E,T] = annihilation2D;

if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf,
  pdegplot(G,'subdomainlabels','on'), axis equal

  figure(2), clf,
  pdemesh(P,E,T), axis equal
end

% simulate
umod = urdme(umod,'seed',20170221);
umod = urdme2pde(umod);

% compute temporal means
mean_A = mean(umod.pde.U(1,:,end/2:end),3);
mean_B = mean(umod.pde.U(2,:,end/2:end),3);

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(3), clf
  pdesurf(umod.pde.P,umod.pde.T,mean_A');
  title('Mean A');
  
  figure(4), clf
  pdesurf(umod.pde.P,umod.pde.T,mean_B');
  title('Mean B')
end
