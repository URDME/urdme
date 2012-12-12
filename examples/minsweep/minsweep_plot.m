function [slices,conc] = minsweep_plot(index)
%MINSWEEP_PLOT Postprocessing script for the minsweep model.
%   MINSWEEP_PLOT(index) Computes and plots the temporal average of MinD_m
%   along with the total number of MinD_m in one half of the
%   bacterium. Results (precomputed) are taken from files #index in
%   results/.
%
%   [SLICES,CONC] = MINSWEEP_PLOT(index), with two output arguments no
%   figures are created. Instead the temporal average CONC in slices SLICES
%   is computed. Hence plot(SLICES,CONC) displays the temporal average of
%   the relative concentration of MinD_m over the bacterium.

% P. Bauer   2012-10-29 (Port to comsol 4.x)
% S. Engblom 2011-11-21

% load umod+u
load(sprintf('results/in%d',index));

if strcmp(umod.mph,'4.x')
  umod.comsol=mphload(sprintf('%s/results/mod%d.mph',pwd,index));
end

%insert solution to comsol model
%umod = urdme2comsol(umod,U,1);
umod = urdme_addsol(umod,sprintf('results/out%d',index));


% compute temporal average
Umean = (umod.U(:,500:end-1)*diff(umod.tspan(500:end))')/ ...
        (umod.tspan(end)-umod.tspan(500));
% (disregard the first 500 steps)

% find x-coordinate and surface/volume of all surface voxels
sdn2 = []; %find(fem.urdme.sd ~= 2); % remove non-boundary elements
% (does note make a difference since MinD_m only lives on the surface)
Mspecies = size(umod.N,1);

if strcmp(umod.mph,'4.x')
  xmi = mphxmeshinfo(umod.comsol);
  x = xmi.dofs.coords(1,1:Mspecies:end); 
else
  if ~isfield(umod.comsol,'xmesh')
    umod.comsol.xmesh = meshextend(umod.comsol);
  end
  try
  dofs = xmeshinfo(umod.comsol,'out','dofs');
  catch
    umod.comsol.xmesh = meshextend(umod.comsol);
    dofs = xmeshinfo(umod.comsol,'out','dofs');
  end
  x = dofs.coords(1,1:Mspecies:end);  
end

xx = x;
xx(sdn2) = [];

% surface/volume
if strcmp(umod.mph,'4.x')
  DM = mphmatrix(umod.comsol,'sol1','out','D');
  M = DM.D;
  vol = full(sum(M,2));  
else
  M = assemble(umod.comsol,'out',{'D'});% surface concentration
  vol = sum(M);
end
vol = vol(1:Mspecies:end);
% uncomment to get volume concentration instead:
%vol = fem.urdme.vol;
vol(sdn2) = [];

% insert first and mean solution into fem2 for visualization
umod2=umod;
% 'fool' Comsol to produce surface concentration:
umod2.vol = vol;

umod2.tspan=umod.tspan(end);
%umod2 = urdme2comsol(umod2,Umean,1); %TODO: why is this here?

% slices in the x-direction
minx = min(x);
maxx = max(x);
slices = linspace(minx,maxx,30);

% count MinD_m (temporal average) and total volume in each slice
ixMinD_m = 2; % index of MinD_m according to URDME's ordering
MinD_m = Umean(ixMinD_m:Mspecies:end);
MinD_m(sdn2) = [];
[~,ix] = histc(xx,slices);
aver = full(sparse(ix,1,MinD_m,size(slices,2),1));
vol = full(sparse(ix,1,vol,size(slices,2),1));
conc = (aver./vol)'; % concentration per slice
conc(isnan(conc)) = 0; % (avoid 0/0)
conc = conc./max(conc); % normalization

% count MinD_m close to the poles
mx = 0.5*(maxx+minx);
lx = maxx-minx;
MinD_m = umod.U(ixMinD_m:Mspecies:end,:);
tMinD_m1 = sum(MinD_m(x <= mx-0.25*lx,:));
tMinD_m2 = sum(MinD_m(x > mx+0.25*lx,:));

%% final plots
%if nargout == 0
%  figure(1),
%  subplot(2,1,1);
%  plot(slices,conc,'b.-','LineWidth',2);
%  set(gca,'YLim',[0 1],'YTick',0.0:0.2:1.0, ...
%          'XLim',[minx maxx],'XTick',[minx 0 maxx]);
%  title('Relative surface concentration (temporal average)');
%  xlabel('x [\mu m]');
%  
%  subplot(2,1,2);
%  plot(umod.tspan,tMinD_m1,'b',umod.tspan,tMinD_m2,'r','Linewidth',2);
%  set(gca,'YLim',[0 1000],'YTick',0:250:1000, ...
%          'XLim',[1000 1500],'XTick',1000:100:1500);
%  title('Copy number near the poles');
%  xlabel('t [s]');
%  ylabel('# molecules');
%  
%  % PB: This doesnt seem to work!
%  figure(2);
%  if strcmp(umod.mph,'4.x')
%    mphplot(umod2.comsol,'Tetdata','MinD_m');
%  else
%    postplot(umod2.comsol,'Tetdata','MinD_m');
%  end
%  title('Surface concentration: temporal average');
%end

% final plots
if nargout == 0
  figure(1),
  subplot(2,1,1);
  plot(slices,conc,'b.-','LineWidth',2);
  set(gca,'YLim',[0 1],'YTick',0.0:0.2:1.0, ...
          'XLim',[minx maxx],'XTick',[minx 0 maxx]);
  ht=title('Average concentration of MinD_m');
  %title('Relative surface concentration (temporal average)');
  hx=xlabel('         location (\mu m)');
  set(hx,'FontSize',12)
  set(ht,'FontSize',16)

  subplot(2,1,2);
  plot(umod.tspan,tMinD_m1,'b',umod.tspan,tMinD_m2,'r','Linewidth',2);
  set(gca,'YLim',[0 1000],'YTick',0:250:1000, ...
          'XLim',[1000 1500],'XTick',1000:100:1500);
  %ht=title('Copy number near the poles');
  ht=title('MinD_m polar oscillations        ');
  hx=xlabel('time (s)');
  hy=ylabel('# molecules');
  set(hx,'FontSize',14)
  set(hy,'FontSize',12)
  set(ht,'FontSize',16)
  legend('Left','Right');


  figure(2); 
  %postplot(fem2,'Tetdata','MinD_m','Geom','off','Tetbar','off');
  if strcmp(umod.mph,'4.x')
    mphplot(umod2.comsol,'Tetdata','MinD_m','Geom','off','Tetbar','off');
  else
    postplot(umod2.comsol,'Tetdata','MinD_m','Geom','off','Tetbar','off');
  end
  %title('Surface concentration: temporal average');
  view(0,0);
  set(gca,'ZTickLabel',[])
  set(gca,'ZTick',[])
  set(gca,'XTickLabel',[])
  set(gca,'XTick',[])
end
