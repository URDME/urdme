% Annihilation example.
%   Species A and B start at opposite boundaries (in 1D) and
%   annihilate when they meet. This example shows how inline
%   propensities are used, and how one can set up a simple 1D
%   simulation.

% S. Engblom 2017-02-18 (Major revision, URDME 1.3, Comsol 5)
% B. Drawert 2012

% generate URDME model
clear umod
umod = annihilation;

% simulate
umod = urdme(umod,'seed',20170219);

% compute means
mean_A = mean(umod.U(1:2:end,end/2:end),2)./umod.vol;
mean_B = mean(umod.U(2:2:end,end/2:end),2)./umod.vol;

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf
  plot(umod.private.mesh,mean_A,'-sb',umod.private.mesh,mean_B,'-dr');
  legend('A','B');
  xlabel('x');
end
