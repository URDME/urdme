%Morphogenesis 3D script file.
%   This file runs the Schnakenberg model in a 3D torus.

% S. Engblom 2017-03-17

% load Comsol model

%% (1) geometry: load Comsol model
if exist('mli','dir')
  model = mphload('torus.mph');

  % create URDME struct
  umod = comsol2urdme(model);
else
  % alternative: simply load the equivalence of the above construct
  load torus
end

% add reactions
umod = schnakenberg(umod);
if exist('tspan','var')
  umod.tspan = tspan;
end

% volume scaling
umod.vol = 50/mean(umod.vol)*umod.vol;

% solve
umod = urdme(umod,'propensities','schnakenberg','report',0);

return;

% visualize using Comsol
umod = urdme2comsol(umod);

figure(1), clf
umod.comsol.result.create('res1','PlotGroup3D');
umod.comsol.result('res1').feature.create('surf1','Surface');
umod.comsol.result('res1').feature('surf1').set('expr','U');

M = struct('cdata',{},'colormap',{});
el = linspace(30,60,numel(umod.tspan));
for i = 1:numel(umod.tspan)
  umod.comsol.result('res1').set('t',num2str(umod.tspan(i)));
  mphplot(umod.comsol,'res1');
  set(gca,'Xtick',[],'Ytick',[],'Ztick',[],'title',[]);
  view(el(i)-7.5,el(i));
  drawnow;
  M(i) = getframe;
end

% create animation
movie2gif(M,{M([1 2 10 20]).cdata},'schnakenberg.gif', ...
          'delaytime',0.1,'loopcount',0);
