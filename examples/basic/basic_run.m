%Basic example script.
%   Reaction X+Y <--> Z in a 3D sphere.

% S. Engblom 2019-11-27 (Revision, augment model using rparse)
% S. Engblom 2017-05-09 (Revision, rparse)
% S. Engblom 2017-02-19 (Revision)

%% (1) geometry: load Comsol model
if exist('mli','dir')
  model = mphload('basic.mph');

  % create URDME struct
  umod = comsol2urdme(model);
else
  % alternative: simply load the equivalence of the above construct
  load sphere
end
  
%% (2) reactions
umod = rparse(umod, ...
    {'X+Y > X*Y/vol > Z' ...
     'Z > Z > X+Y'}, ...
    {'X' 'Y' 'Z'},{},'basic.c');

% initial number of species
Mspecies = 3; % ordering is [X Y Z]
Ncells = numel(umod.vol); 
umod.u0 = zeros(Mspecies,Ncells); 

% distribute 50 X and 50 Y randomly
rng(123); % to get reproducible results
cell = randi(Ncells,1,50);
umod.u0(1,:) = full(sparse(1,cell,1,1,Ncells));
cell = randi(Ncells,1,50);
umod.u0(2,:) = full(sparse(1,cell,1,1,Ncells));

%% (3) simulate

% diffusion is set to 0.1
umod = urdme(umod,'tspan',0:.1:100,'seed',123);

%% (4) postprocessing using Comsol
if exist('mli','dir') && (~exist('plotting_off','var') || ~plotting_off)
  umod = urdme2comsol(umod);  % put result back into Comsol

  figure(1), clf, mphmesh(model);

  figure(2), clf,
  umod.comsol.result.create('res1','PlotGroup3D');
  umod.comsol.result('res1').set('t','1');
  umod.comsol.result('res1').feature.create('surf1','Surface');
  umod.comsol.result('res1').feature('surf1').set('expr','Z');
  mphplot(umod.comsol,'res1');

  figure(3), clf,
  umod.comsol.result.create('res2','PlotGroup3D');
  umod.comsol.result('res2').set('t',sprintf('%d',umod.tspan(end)));
  umod.comsol.result('res2').feature.create('surf1','Surface');
  umod.comsol.result('res2').feature('surf1').set('expr','Z');
  mphplot(umod.comsol,'res2');
end

% extract Z from the result
Z_series_fastDiff = umod.U(3:3:end,:);

%% (5) re-iterate and run again

if exist('mli','dir') 
  % change diffusion coefficients to .001
  model.physics('chds').feature('cdm1').set('D_X',{'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});
  model.physics('chds').feature('cdm1').set('D_Y',{'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});
  model.physics('chds').feature('cdm1').set('D_Z',{'.001[m^2/s]' '0' '0' '0' '.001[m^2/s]' '0' '0' '0' '.001[m^2/s]'});

  % run the modified example
  vmod = comsol2urdme(model);
  % keep only these fields:
  umod.D = vmod.D;
  umod.comsol = vmod.comsol;
  % (the rest of umod is re-used)
else
  % simple way:
  umod.D = 0.01*umod.D;
end
umod = urdme(umod,'compile',0); % already compiled

%% (6) postprocessing using Comsol
if exist('mli','dir') && (~exist('plotting_off','var') || ~plotting_off)
  umod = urdme2comsol(umod);

  figure(4), clf,
  mphplot(umod.comsol,'res1');

  figure(5), clf,
  mphplot(umod.comsol,'res2');
end

% extract Z
Z_series_slowDiff = umod.U(3:3:end,:);

%% (7) visualize
if ~exist('plotting_off','var') || ~plotting_off
  figure(7), clf,
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
