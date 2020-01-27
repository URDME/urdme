% The MinCDE model solved with cytosol subdiffusion.

% S. Engblom 2017-04-12 (Major revision)
% S. Engblom 2017-02-20 (Revision)
% S. Engblom 2014-05-14

%% (1) geometry: load Comsol model
if exist('mli','dir')
  model = mphload('../../examples/mincde/coli.mph');

  % create URDME struct
  umod = comsol2urdme(model);
  
  % compute h and save it at ldata
  xmi = mphxmeshinfo(umod.comsol);
  Mspecies = numel(xmi.fieldnames);
  umod.ldata = mphinterp(umod.comsol,'h','coord', ...
                         xmi.dofs.coords(:,1:Mspecies:end),'solnum',1);
  umod.ldata = umod.ldata(xmi.dofs.nodes(1:Mspecies:end)+1);

  % also save the x-coordinates (for postprocessing)
  umod.private.x = xmi.dofs.coords(1,1:Mspecies:end);
else
  % alternative: simply load the equivalence of the above construct
  load ../../examples/mincde/coli
end

% participating species
species = {'MinDcytATP' 'MinDmem' 'MinE' 'MinDE' 'MinDcytADP'};
Mspecies = numel(species);
% Species MinDmem (#2) and MinDE (#4) only diffuse on the membrane,
% while species MinDcytATP (#1), MinE (#3), and MinDcytADP (#5) are
% subject to subdiffusive (volume) crowding.

%% (2) form the (sub-)diffusion operator

% remove the elements corresponding to the cytosol for the species
% bound to the membrane
D = umod.D;
cyt = find(umod.sd == 1);
ixremove = zeros(0,1);
for s = [2 4]
  ixremove = [ixremove; Mspecies*(cyt-1)+s];
end
D(ixremove,:) = 0;
D(:,ixremove) = 0;
d = full(sum(D,1));
D = D+sparse(1:size(D,1),1:size(D,1),-d);

% fetch internal states approximation from pre-simulations
[gamma,A,p] = subdata(0.2);
Nstates = numel(gamma);

% initial- and equilibrium effective diffusion
[gamma0,igammamax] = max(gamma); % fastest diffusing rate/state
gammabar = p'*gamma;             % equilibrium diffusion

% note: scaling with gammabar to agree at steady-state
select = logical([1 0 1 0 1]);
umod.D = subdmatrix(A,gamma/gammabar,D,select);
state1 = cumsum([1 1+select*(Nstates-1)]);
% the 5 species are now found in [1..10, 11, 12..21, 22, 23..32], and
% state1(i) points to the first state of species i

% initial data: sample a certain number of species 3 and 5, roughly
% uniformly spatially distributed and all in the fastest diffusing
% state
Ncells = numel(umod.sd);
u0 = zeros(state1(end)-1,Ncells);
nMinE = 1040;
nMinD = 4002;

ind = floor(Ncells*rand(1,nMinE))+1;
u0(state1(3)+igammamax-1,:) = full(sparse(1,ind,1,1,Ncells));

ind = floor(Ncells*rand(1,nMinD))+1;
u0(state1(5)+igammamax-1,:) = full(sparse(1,ind,1,1,Ncells));
umod.u0 = u0;

%% (3) form the subdiffusive reactions

% sequence expansion: expand all species in i = 1:Nstates inner states
seq = {['ie']' [1:Nstates; 0:Nstates-1]};
% (e = 0:Nstates-1 is used to expand kdD[e] due to C-arrays)

% append '_$i' to species names with internal states
species(select) = cellfun(@(x)cat(2,x,'_$i'), ...
                          species(select),'uniformoutput',0);

% reactions, with internal states marked out
r1 = 'MinDcytATP_$i > sd == 2 ? kd*MinDcytATP_$i/ldata[0] : 0.0 > MinDmem';
r2 = ['MinDcytATP_$i + MinDmem > ' ...
      'kdD[$e]*MinDcytATP_$i*MinDmem/(1000.0*NA*vol) > ' ...
      'MinDmem+MinDmem'];
r3 = 'MinE_$i + MinDmem > kde*MinE_$i*MinDmem/(1000.0*NA*vol) > MinDE';
r4 = 'MinDE > ke*MinDE > MinDcytADP_10 + MinE_10';
r5 = 'MinDcytADP_$i > kp*MinDcytADP_$i > MinDcytATP_10';
% ('10' is igammamax, written out explicitly for clarity)

% reaction rates
kdD = 9.0e6*1.65;
% (the *1.65 takes us to a slightly more unstable regime compared to
% the original)
kde = 5.56e7;
NA = 6.022e23;
% varying but normalized rates:
weights = linspace(1,1.3,Nstates);
weights = weights/(weights*p);
rates = {'NA' 6.022e23 'kd' 1.25e-8 'ke' 0.7 'kp' 0.5 ...
         'kdD' kdD*weights 'kde' kde};
[~,umod.N,umod.G] = rparse({r1 r2 r3 r4 r5},species,rates, ...
                           'mincde_subdiff.c',seq);

% output times
umod.tspan = 0:200;

% solve
umod.seed = 321;
umod = urdme(umod,'propensities','mincde_subdiff.c','report',2);

return;

% visualization using pre-stored simulations

% original run, but with kdD =  9.0e6*1.65
load data/originalrun

% fetch the x-coordinates
xx = umod.private.x;
minx = min(xx);
maxx = max(xx);

% count MinDmem close to the poles over time
mx = 0.5*(maxx+minx);
lx = maxx-minx;
% it is interesting to observe also the transient here:
Mspecies = size(umod.N,1);
ixMinDmem = 2;
MinDmem = umod.U(ixMinDmem:Mspecies:end,:);
tMinDm1 = sum(MinDmem(xx <= mx-0.25*lx,:));
tMinDm2 = sum(MinDmem(xx > mx+0.25*lx,:));

figure(1), clf
subplot(2,1,1);
plot(umod.tspan,tMinDm1,'b',umod.tspan,tMinDm2,'r--','Linewidth',2);
%title('MinDmem polar oscillations');
ylabel('# molecules');
axis tight
ax = axis;

figure(2), clf
F = fftshift(fft(tMinDm1-mean(tMinDm1)));
f = (-100:100)/numel(umod.tspan);
F = abs(F);
semilogy(f(2:5:end),F(2:5:end),'b','Linewidth',2);
%title('Power spectrum');
axis([0 f(end) 1e2 1e5])
ylabel('Power'); xlabel('f');
hold on

% example using subdiffusion (the above code)
load data/subdrun

Nstates = 10;
Ncells = numel(umod.sd);
select = logical([1 0 1 0 1]);
state1 = cumsum([1 1+select*(Nstates-1)]);
U = reshape(umod.U,state1(end)-1,Ncells,[]);

% x-coordinates
xx = umod.private.x;
minx = min(xx);
maxx = max(xx);

% count MinDmem
ixMinDmem = 2; % index of MinDmem according to URDME's ordering
MinDmem = squeeze(sum(U(state1(ixMinDmem):state1(ixMinDmem+1)-1,:,:),1));

% count MinDmem close to the poles over time
mx = 0.5*(maxx+minx);
lx = maxx-minx;
tMinDm1 = sum(MinDmem(xx <= mx-0.25*lx,:));
tMinDm2 = sum(MinDmem(xx > mx+0.25*lx,:));

figure(1),
subplot(2,1,2);
plot(umod.tspan,tMinDm1,'b',umod.tspan,tMinDm2,'r--','Linewidth',2);
xlabel('time (s)');
ylabel('# molecules');
axis(ax);

figure(2), 
F = fftshift(fft(tMinDm1-mean(tMinDm1)));
f = (-100:100)/numel(umod.tspan);
F = abs(F);
semilogy(f(2:5:end),F(2:5:end),'b--','Linewidth',2);
axis([0 f(end) 1e2 1e5])
ylabel('Power'); xlabel('f');

return;

% possibly print to files
figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
print -depsc ~/oscillations.eps

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 250]);
print -depsc ~/power.eps
