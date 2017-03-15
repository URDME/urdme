%Basic example script.
%   Reaction X+Y <--> Z in a 3D sphere.

% S. Engblom 2017-02-19 (Revision)

% load Comsol file
model = mphload('basic.mph');

% display mesh
figure(1), clf,
mphmesh(model);

% create URDME struct
umod = comsol2urdme(model); % get geometry and diffusion from Comsol
umod = basic(umod);         % add reactions

% run this example (diffusion is set to 0.1)
umod = urdme(umod,'propensities','basic');

umod = urdme2comsol(umod);  % put result back into Comsol

% display result using Comsol
figure(2), clf,
umod.comsol.result.create('res1','PlotGroup3D');
umod.comsol.result('res1').set('t','1');
umod.comsol.result('res1').feature.create('surf1', 'Surface');
umod.comsol.result('res1').feature('surf1').set('expr', 'Z');
mphplot(umod.comsol,'res1');

figure(3), clf,
umod.comsol.result.create('res2','PlotGroup3D');
umod.comsol.result('res2').set('t',sprintf('%d',umod.tspan(end)));
umod.comsol.result('res2').feature.create('surf1', 'Surface');
umod.comsol.result('res2').feature('surf1').set('expr', 'Z');
mphplot(umod.comsol,'res2');

% extract Z from the result
Z_series_fastDiff = umod.U(3:3:end,:);

% change diffusion coefficients to .001
model.physics('chds').feature('cdm1').set('D_X', {'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});
model.physics('chds').feature('cdm1').set('D_Y', {'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});
model.physics('chds').feature('cdm1').set('D_Z', {'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});

% run the modified example
umod = comsol2urdme(model);
umod = basic(umod);
umod = urdme(umod,'compile',0); % already compiled

% plot
umod = urdme2comsol(umod);

figure(4), clf,
mphplot(umod.comsol,'res1');

figure(5), clf,
mphplot(umod.comsol,'res2');

% extract Z
Z_series_slowDiff = umod.U(3:3:end,:);

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(6), clf,
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
