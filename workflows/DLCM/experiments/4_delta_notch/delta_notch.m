%DELTA_NOTCH Delta Notch signalling in a growing population.
%   This model runs the cell population physics from BASIC_TEST but
%   with a simple proliferation process and a delta-notch model at the
%   same time (both in continuous time). The purpose of the experiment
%   is to show how continuous-in-time processes may be easily coupled
%   together.
%
%   Handles both discrete (default) or continous quantities.
%   ~2-3 min simulation time on a laptop.

% E. Blom 2024-11-12 (major revision, using URDME's DLCM solver)
% S. Engblom 2017-12-28 (minor revision)
% S. Engblom 2017-08-30 (minor revision)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-27 (finalized physics)
% S. Engblom 2016-12-21 (updated the physics, changed to hex-cells)
% S. Engblom 2016-11-10

%% (1) geometry
Nvoxels = 41;
mesh_type = 2;  % cartesian mesh

% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 5000;
end
Tres = 100;
if ~exist('seed', 'var')
  seed = 123;
end
if ~exist('internal_state', 'var')
  internal_state = 'discr';   % {'discr'} | 'cont'
end
ntypes = 1; % number of cell types: 1: living cells

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.

% get boundary dofs, extdof
hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
extdof = find(P(1,:).^2 + P(2,:).^2 > (0.866 - ...
hmax)^2 & P(1,:).^2 + P(2,:).^2 <= 0.866^2);

%% (2) State transitions per cell

% Input leged: Q is an n-by-Nvoxels matrix, where n is the number of
% Laplace quantities. Q(:,i=1) is reserved for pressure, and generally
% for i = 1..n. Inputs are U = cell number, Q = scalar field,
% P = dof2position map, T = time. QI is an Nvoxels-by-2-by-
% -Ninternal local cell data array, where 2 is maximum number of
% cells per voxel and Ninternal is number of internal states (voxel index
% and cell type as 1st and 2nd position, respectively).

% Rates is a function that returns a 1-by-nmig cell array, where n is the
% number of migration 'potentials': migration rates are proportional to the
% potential gradients. Inputs to are U = cell number, Q = scalar field, QI
% is Ncells-by-Ncapacity-by-ninternal local cell data array,
% P = dof2position map, T = time.
Rates = @(U,Q,QI,P,t){Q(:,1)};  % here, migration

% Drate is a function that returns a 1-by-nmig cell array that
% specifies how the corresponding migration rates are scaled. Uf is the
% cell number in the voxel to move from, Ut in voxel to move to.
% Note: solver considers only migration from bdof_m and sdof_m!
Drate = @(Uf,Ut,Q,QI,P,t){1.*(Uf==1).*(Ut==0)+1.*(Uf==2).*(Ut<2)};

%% (3) Form population

% initial small population
r0 = 0.2;
ii1 =  find(abs(P(1,:)) < r0 & abs(P(2,:)) <= r0);   % alive cells
% U is ntype-by-ncells sparse vector, representing the cell population
U = fsparse(ii1(:),1,2,[Nvoxels^2 1]); % doubly occupied

%% (4) "outer" URDME-struct instead:
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% 'UL' means cell of phenotype L.
% Define all reaction events first, and after that quantities
mu_prol = 0.01;                                      % proliferation rate
umod = rparse(umod, ...
              {'U1 > mu_prol*(U1 == 1) > U1+U1', ... % proliferation event
              'Q1 > (U1>1) > Q1+Q1'}, ... % pressure sources
              {'U1' 'Q1'}, ...
              {'mu_prol' mu_prol}, ...
              'delta_notch_outer');
umod.u0 = [full(U)'; zeros(1,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps

%% (5) Delta-Notch model with random IC
% Delta-Notch (D-N) model parameters
a = 0.01;
b = 100;
v = 1;
k = 2;
h = 2;
% delta-notch timescale: 0.01 is "slow" compared to mechanics, 1 or
% above is "fast"
Tdn = 0.04;

% random (0,Dmax or Nmax) initial delta & notch values in each cell:
rng(123);
nstates = 2;                 % nr of distinct internal states
Dmax = 400; Nmax = Dmax;
ii2 = ii1;                   % (every voxel has 2 cells)
De = fsparse(ii1(:),1,randi(Dmax,numel(ii1),1),[Nvoxels^2 2]);
De(ii2,2) = randi(Dmax,numel(ii2),1);
No = fsparse(ii1(:),1,randi(Nmax,numel(ii1),1),[Nvoxels^2 2]);
No(ii2,2) = randi(Nmax,numel(ii2),1);

% celldata is an Ntypes*2-by-Ncells array of cell state data -> mumod.u0.
% Here containing Delta and Notch values for each cell
celldata = zeros(nstates*2, Nvoxels^2);
adof = find(sum(U,2)>0); sdof = find(sum(U,2)>1);
celldata(1, adof) = De(adof,1);   % first cell in voxel
celldata(2, sdof) = De(sdof,2);   % second cell in voxel
celldata(3, adof) = No(adof,1);
celldata(4, sdof) = No(sdof,2);

%% (6) Internal States using URDME SSA
% define reactions
r1 = '@ > (sd >= $i)*Tdn*v/(1+b*pow(N$i/Nmax,h))*Nmax > D$i';
r2 = 'D$i > Tdn*D$i*v > @';
r3 = '@ > (sd >= $i)*Tdn*ldata[0]*Nmax > N$i';
r4 = 'N$i > Tdn*N$i > @';

% Maximum micro-time step (during which ldata is fixed) by "5%-norm rule"
maxdt_fun = @(U, Q, QI, ldata_fun)0.05*norm([QI(:,:,3); QI(:,:,4)],'fro')...
  /(Tdn*norm([v./(1+b*(QI(:,:,4)/Nmax).^h)*Nmax;-QI(:,:,3)*v; ...
  ldata_fun{1}*ones(1,2);-QI(:,:,4)],'fro'));

% Define reaction rates ldata[1:2] (held constant during maxdt)
F = @(x)(x.^k./(a+x.^k));         % (simplifies expression below)
% Function for ldata input that returns a 1-by-nldata cell array, with
% ldata_fun(...){l} containing ldata[l] for l = 1:nldata.
% Here the function return the avg. neighboring Delta activity (ignoring
% cells in same voxel; can be accounted for by altering r3)
ldata_fun = @(U, Q, QI, Ne){F(sum(Ne*QI(:,:,3),2)/Nmax./max(Ne*U,1))};

% Construct, parse, and compile ssa solver fully outside dlcm
Dexpr = cell(2*nstates,1); Dexpr(1:2*nstates) = {1};
mumod = pde2urdme(P,T,Dexpr);
mumod = rparse(mumod,{r1 r2 r3 r4},{'D$i' 'N$i'}, ...
{'v', v, 'b' b 'h' h 'Tdn' Tdn 'Nmax', Nmax}, ...
'delta_notch_inner', {'i', 1:2});
mumod.seed = seed;
Nlive = size(De,1);                     % nr of live cells
Nspec = nstates*2;
mumod.D = sparse(Nspec*Nlive,Nspec*Nlive);
mumod.u0(1:Nspec,:) = celldata;         % internal states
mumod.vol = ones(1, Nlive);
mumod.sd = ones(1, Nlive);              % compile dummy
mumod.tspan = [0 1];                    % compile dummy
mumod = urdme(mumod, 'solve', 0);       % parse & compile

% finally load essentials and mumod with internal specifics into umod
umod = dlcm2urdme(umod, P, gradquotient, [],[],[], 'Rates', Rates, ...
  'Drate', Drate, 'mumod', mumod, 'ldata_fun', ldata_fun, ...
  'maxdt_fun', maxdt_fun, 'internal_state', internal_state);

%% (7) solve
umod = urdme(umod,'solver','dlcm', 'solve',1, 'seed', 123);
