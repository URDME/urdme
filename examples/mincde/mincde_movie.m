%MINCDE_MOVIE Create animation of Min-model.

% S. Engblom 2017-03-15

% assumes mincde_run was successful
umod = urdme2comsol(umod);

figure(3), clf
umod.comsol.result.create('res2','PlotGroup3D');
umod.comsol.result('res2').feature.create('surf1','Surface');
umod.comsol.result('res2').feature('surf1').set('expr','MinDmem');

M = struct('cdata',{},'colormap',{});
el = linspace(-20,20,numel(umod.tspan));
for i = 1:numel(umod.tspan)
  umod.comsol.result('res2').set('t',num2str(umod.tspan(i)));
  mphplot(umod.comsol,'res2');
  set(gca,'Xtick',[],'Ytick',[],'Ztick',[],'title',[]);
  view(2*el(i)+3.5,el(i));
  drawnow;
  M(i) = getframe;
end

return;

% create animation
movie2gif(M,{M([1 2 10 20]).cdata},'mincde.gif', ...
          'delaytime',0.1,'loopcount',0);
