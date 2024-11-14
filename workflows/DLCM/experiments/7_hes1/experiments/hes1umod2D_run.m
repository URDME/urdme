%HES1UMOD2D_RUN Runs the URDME Hes1 model in 2D.
%
%   This script also helps building models:
%
%   - Define solve = 0 to build a model without solving.
%
%   - Define VOL to be the cell volume (default: 50).
%
%   - Define Nvoxels == '3' (note: yes, a string) to build
%   specifically for a 3-cell setup.
%
%   See also HES1UMOD.

% S. Engblom 2024-04-10

% turn off to avoid solving and using the script as a model builder
if ~exist('solve','var')
  solve = 1;
end

% turn on to save data
if ~exist('save_data','var')
  save_data = 0;
end

% critical parameter: volume of one voxel
if ~exist('VOL','var')
  % cell volume of mouse embryonal stem cell ~50 um^3
  VOL = 50; % unit: um^3, so "50" is ~one cell 
end
VOL % echo it

% build the geometry and model
clear umod

if ~exist('Nvoxels','var')
  % cells live in a rectangular Nvoxels-by-Nvoxels region
  Nvoxels = 20; % (edit default here)
end
if ~strcmp(Nvoxels,'3')
  % fetch discretization
  [P,E,T] = basic_mesh(2,Nvoxels);
  [V,R] = mesh2dual(P,E,T,'voronoi');

  % (boolean) neighbor matrix
  N_op = dt_neighe(V,R);
  N_op = N_op ~= 0;
else
  % special case: 3-cell problem
  Nvoxels = 3; % note: you must intentionally set Nvoxels = '3'
  N_op = [0 1 1; 1 0 1; 1 1 0];
end
Nvoxels % echo it

% find layers in mesh (for visualization)
layers = mesh_layers(N_op,floor(Nvoxels/2+1));

% use N_op to form a diffusion operator D
D = (N_op./sum(N_op,2))';
D = D-speye(size(N_op,1));

% scaling and cell/voxel volume
rng(20241021);
% note: concentrations with an estimated spread
conc = hes1_conc(size(N_op,1)); % given in micromolar (uM)
avogadro = 6.022*1e23;
% 1 VOL = 1 um^3 - scaling to get volume in litres
cell_vol = VOL*1e-15;
umod.vol = repmat(cell_vol,1,size(N_op,1));
umod.sd = ones(1,size(N_op,1));

% fetch the reactions
umod = hes1umod(umod);
umod.solve = solve;

% go from initially given concentrations to # molecules:
%   moles [mol] = molar concentration [M] * volume [L]
%   # molecules = moles [mol] * avogadro's number [#/mol]
scaling_factor = conc*1e-6*cell_vol*avogadro;

% initial values (close to given concentrations but in # molecules)
u0 = umod.u0.*scaling_factor([1 1:5],:); % D = Din
% "stochastic rounding":
%umod.u0 = floor(u0)+(rand(size(u0)) < u0-floor(u0));
% rounding:
umod.u0 = round(u0);

% define the diffusion matrix: the rate alphaN is interpreted as a
% diffusion event - albeit a special one
alphaN = umod.private.rp.RateVals{strcmp(umod.private.rp.RateNames,'alphaN')};
% incorporate also scaling by SCALE:
SCALE = umod.gdata;
D = alphaN/SCALE*D;

% the negative (outgoing) part of D concerns species Din, the positive
% (incoming) part of D concerns N, i.e., Din diffuses into a
% neighboring voxel, transforming into an N at the same time
ixDin = strcmp(umod.private.rp.Species,'Din');
ixN = strcmp(umod.private.rp.Species,'N');
% perhaps not 100% obvious, but this is how it is done:
Mspecies = size(umod.N,1);
umod.D = kron(D.*(D < 0),diag(sparse(ixDin)))+ ...
         kron(D.*(D > 0),sparse(find(ixN),find(ixDin),1,Mspecies,Mspecies));

% solve
disp('Solving with NSM...');
umod = urdme(umod,'solver','nsm','report',2);
disp('   ...done.');

disp('Solving with UDS...');
vmod = umod;
vmod.solverargs = ['odesolv' {@ode15s} ...
                   'odeopts' odeset('RelTol',1e-4,'AbsTol',1e-2)];
vmod = urdme(vmod,'solver','uds');
disp('   ...done.');

if save_data
  disp('Saving...')
  save("../data/hes1umod_vol"+num2str(VOL)+".mat",'V','R','umod','VOL', ...
       'layers','N_op')
  disp('...saved.')
end

return;

% plot mean in layers and per fate
fig = 0;
for wmod = [umod vmod]
  fig = fig+1;
  % change from # molecules back to concentrations in uM
  P = wmod.U(5:Mspecies:end,:)./wmod.vol'/avogadro*1e6;
  tspan = wmod.tspan;
  
  % judge end-fate
  [Pincr,ix] = sort(P(:,end));
  [~,ijmp] = max(diff(Pincr));
  Plo = nan(size(P));
  Plo(ix(1:ijmp),:) = P(ix(1:ijmp),:);
  Phi = nan(size(P));
  Phi(ix(ijmp+1:end),:) = P(ix(ijmp+1:end),:);

  % final plot
  figure(fig), clf,
  cmap = colormap(hsv(floor(Nvoxels/2+1)));
  for i = 1:size(layers,1)
    plot(tspan/60,nanmean(Plo(layers{i,2},:)),'Color',cmap(i,:)), 
    hold on,
    plot(tspan/60,nanmean(Phi(layers{i,2},:)),'Color',cmap(i,:)), 
  end
  xlabel('time [h]')
  ylabel('P in all cells [um]');
  axis([0 tspan(end)/60 0 1]);
end
