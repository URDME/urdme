%Script for the MinCDE model.
%   Setups and simulates the Min system. Postprocessing in the form of
%   temporal average of MinDmem along with the total number of MinDmem
%   near the poles.

% S. Engblom 2017-02-20 (Revision, URDME 1.3, Comsol 5.2)
% P. Bauer 2012-09-03
% A. Hellander 2010-06-09 (temporal average)
% J. Cullhed 2008-08-14

% load Comsol model
model = mphload('coli.mph');

% create URDME struct
umod = comsol2urdme(model);
umod = mincde(umod); % add reactions
if exist('tspan','var')
  umod.tspan = tspan;
end

% run model
umod = urdme(umod,'propensities','fange','report',0);

% postprocessing, temporal average
U = umod.U;
U2 = mean(U(:,floor(end/2):end),2); % disregard the transient phase

% use a copy of umod to plot
umod2 = umod;
umod2.tspan = umod2.tspan(end);
umod2.U = U2;
umod2 = urdme2comsol(umod2);

if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf
  umod2.comsol.result.create('res1','PlotGroup3D');
  umod2.comsol.result('res1').feature.create('surf1','Surface');
  umod2.comsol.result('res1').feature('surf1').set('expr','MinDmem');
  mphplot(umod2.comsol,'res1');
  title('Temporal average');
  set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
end

% slice averages: fetch some information about the mesh
Mspecies = size(umod.N,1);
xmi = mphxmeshinfo(umod.comsol);
xx = xmi.dofs.coords(1,1:Mspecies:end); 
minx = min(xx);
maxx = max(xx);
slices = linspace(minx,maxx,40);

% count MinDmem (temporal average) and total volume in each slice
ixMinDmem = 2; % index of MinDmem according to URDME's ordering
MinDmem = U2(ixMinDmem:Mspecies:end);

% determine to what slice each voxel belongs to
[~,ix] = histc(xx,slices);

% sum up voxel volumes and molecular counts and determine the temporal
% average
aver = full(sparse(ix,1,MinDmem,size(slices,2),1));
vol = full(sparse(ix,1,umod.vol,size(slices,2),1));
conc = (aver./vol)'; % concentration per slice
conc = conc./max(conc); % normalization

% count MinDmem close to the poles over time
mx = 0.5*(maxx+minx);
lx = maxx-minx;
% it is interesting to observe also the transient here:
MinDmem = umod.U(ixMinDmem:Mspecies:end,:);
tMinDm1 = sum(MinDmem(xx <= mx-0.25*lx,:));
tMinDm2 = sum(MinDmem(xx > mx+0.25*lx,:));

if ~exist('plotting_off','var') || ~plotting_off
  figure(2), clf
  subplot(2,1,1);
  plot(slices,conc,'b.-','LineWidth',2);
  title('Average concentration of MinDmem');
  xlabel('location (\mu m)');
  axis tight

  subplot(2,1,2);
  plot(umod.tspan,tMinDm1,'b',umod.tspan,tMinDm2,'r','Linewidth',2);
  title('MinDmem polar oscillations');
  xlabel('time (s)');
  ylabel('# molecules');
  axis tight
end

return;

% possibly print to files

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
%print -depsc  ../../doc/fig/min_average.eps

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
%print -depsc  ../../doc/fig/min_plots.eps
