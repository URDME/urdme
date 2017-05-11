%Script for the MinCDE model.
%   Setups and simulates the Min system. Postprocessing in the form of
%   temporal average of MinDmem along with the total number of MinDmem
%   near the poles.

% S. Engblom 2017-04-20 (new workflow, offline version)
% S. Engblom 2017-02-20 (Revision, URDME 1.3, Comsol 5.2)
% P. Bauer 2012-09-03
% A. Hellander 2010-06-09 (temporal average)
% J. Cullhed 2008-08-14

%% (1) geometry: load Comsol model
if exist('mli','dir')
  % load Comsol model
  model = mphload('coli.mph');

  % create URDME struct
  umod = comsol2urdme(model);
  
  % One of the propensity functions needs to be scaled by 1/h, where h
  % is the local lengthscale of each subvolume. We can obtain that
  % information by prompting Comsol. Using the ldata vector (input to
  % the propensity functions) we pass h to each propensity function.
  %
  % "mphinterp" is a built-in Comsol function that evaluates an
  % expression at a set of specified points. The predefined expression
  % 'h' can be used to obtain the local length of each subvolume.
  xmi = mphxmeshinfo(umod.comsol);
  Mspecies = numel(xmi.fieldnames);
  umod.ldata = mphinterp(umod.comsol,'h','coord', ...
                         xmi.dofs.coords(:,1:Mspecies:end),'solnum',1);
  umod.ldata = umod.ldata(xmi.dofs.nodes(1:Mspecies:end)+1);

  % also save the x-coordinates (for postprocessing)
  umod.private.x = xmi.dofs.coords(1,1:Mspecies:end);
else
  % alternative: simply load the equivalence of the above construct
  load coli
  Mspecies = 5;
end

%% (2) form the diffusion operator

% We need to modify the diffusion matrix such that the membrane bound
% species only diffuse on the membrane. We achieve this by zeroing out
% the elements corresponding to diffusion in the cytosol.

% use the subdomain vector sd to find the dofs (degress of freedom)
% that are in the cytosol (cyt)
D = umod.D;
cyt = find(umod.sd == 1);

% for MinDmem (#2) and MinDE (#4), remove all dofs in the cytosol
ixremove = zeros(0,1);
for s = [2 4]
  ixremove = [ixremove; Mspecies*(cyt-1)+s];
end
D(ixremove,:) = 0;
D(:,ixremove) = 0;

% enforce sum of rows criterion
d = full(sum(D,1));
D = D+sparse(1:size(D,1),1:size(D,1),-d);
umod.D = D;

% set "random" initial distribution; this can be achieved in several
% ways, below is one example of doing it
Ncells = numel(umod.sd);
u0 = zeros(Mspecies,Ncells);
nMinE = 1040;
nMinD = 4002;

ind = floor(Ncells*rand(1,nMinE))+1;
u0(3,:) = full(sparse(1,ind,1,1,Ncells));

ind = floor(Ncells*rand(1,nMinD))+1;
u0(5,:) = full(sparse(1,ind,1,1,Ncells));

umod.u0 = u0;

% the state of the trajectory will be stored at all time points in
% tspan
if exist('tspan','var')
  umod.tspan = tspan;
else
  umod.tspan = 0:200;  
end

%% (3) form the reactions
r1 = 'MinDcytATP > sd == 2 ? kd*MinDcytATP/ldata[0] : 0.0 > MinDmem';
r2 = ['MinDcytATP + MinDmem > kdD*MinDcytATP*MinDmem/(1000.0*NA*vol) > ' ...
      'MinDmem+MinDmem'];
r3 = 'MinE + MinDmem > kde*MinE*MinDmem/(1000.0*NA*vol) > MinDE';
r4 = 'MinDE > ke*MinDE > MinDcytADP + MinE';
r5 = 'MinDcytADP > kp*MinDcytADP > MinDcytATP';

species = {'MinDcytATP' 'MinDmem' 'MinE' 'MinDE' 'MinDcytADP'};
rates = {'NA' 6.022e23 'kd' 1.25e-8 'kdD' 9.0e6 'kde' 5.56e7 ...
         'ke' 0.7 'kp' 0.5};
[~,umod.N,umod.G] = rparse({r1 r2 r3 r4 r5},species,rates,'fange.c');

% run model
umod = urdme(umod,'propensities','fange','report',0);

%% (4) postprocessing, temporal average
U = umod.U;
U2 = mean(U(:,(end-1)/2:end),2); % disregard the transient phase

% postprocessing using Comsol
if exist('mli','dir') && (~exist('plotting_off','var') || ~plotting_off)
  % use a copy of umod to plot
  umod2 = umod;
  umod2.tspan = umod2.tspan(end);
  umod2.U = U2;
  umod2 = urdme2comsol(umod2);  

  figure(1), clf
  umod2.comsol.result.create('res1','PlotGroup3D');
  umod2.comsol.result('res1').feature.create('surf1','Surface');
  umod2.comsol.result('res1').feature('surf1').set('expr','MinDmem');
  mphplot(umod2.comsol,'res1');
  title('Temporal average');
  set(gca,'Xtick',[],'Ytick',[],'Ztick',[]);
end

% slice averages: fetch the x-coordinates
xx = umod.private.x;
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
