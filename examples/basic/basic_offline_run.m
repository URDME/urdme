%Basic example script. Offline (no Comsol) version.
%   Reaction X+Y <--> Z in a 3D sphere.

% S. Engblom 2017-02-19 (Revision)

clear umod

% load the umod-model where diffusion is set to 0.1
load('offline/offline_d0p01.mat');
umod = basic(umod);

% run this example
umod = urdme(umod,'propensities','basic','seed',29);

% extract Z from the result
Z_series_fastDiff = umod.U(3:3:end,:);

% load the umod-model where diffusion is set to 0.01
load('offline/offline_d0p001.mat')
umod = basic(umod);

% run the modified example
umod = urdme(umod,'compile',0,'seed',29);

% extract Z from the result
Z_series_slowDiff = umod.U(3:3:end,:);

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  subplot(2,1,1)
  plot(sum(Z_series_fastDiff))
  xlim([0 1000])
  ylim([0 50])
  set(gca,'XTickLabel',0:10:100) 
  title('Diffusion = 0.1')
  ylabel('#Z')
  subplot(2,1,2)
  plot(sum(Z_series_slowDiff))
  xlim([0 1000])
  ylim([0 50])
  set(gca,'XTickLabel',0:10:100) 
  title('Diffusion = 0.001')
  ylabel('#Z')
  xlabel('Time')
end
