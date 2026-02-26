%HES1_growth 2D Hes1 model coupled to a proliferation model in DLCM
%   This model runs the cell population physics from BASIC_TEST but
%   with a simple proliferation process and a hes1-notch model at the
%   same time.
%
%   Handles both discrete (default) or continous quantities.
%   cont. <1 min; disc. ~10 min simulation time on a laptop.
%
%   See also HES1UMOD.

% E. Blom 2025-08-14 (coupled to DLCM growth model)
% S. Engblom 2024-04-10

% critical parameter: volume of one voxel
% Set to one for continuous version
if ~exist('VOL','var')
  % cell volume of mouse embryonal stem cell ~50 um^3
  VOL = 50; % unit: um^3, so "50" is ~one cell
end

% build the geometry and model
clear mumod

%% (1) geometry
Nvoxels = 41;
mesh_type = 2;  % hexagonal mesh

% proliferation rate, sets final simulation time
if ~exist('mu_prol', 'var')
  mu_prol = 1/(20*60);  % Expected 20h per cell division
end
% echo it
mu_prol

% Simulate to Tend and save states at Tres intervals
if mu_prol == 0
  R0 = 0.4472;      % initial population radius (static)
else
  R0 = 0.2;         % initial population radius (growing)
end
Tend = 84*60;
Tres = 100;
ntypes = 1; % number of cell types: living cells

if ~exist('seed', 'var')
  seed = 123;
end
if ~exist('internal_state', 'var')
  internal_state = 'discr';   % {'discr'} | 'cont'
end

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - for visualisation

% get boundary dofs, extdof
hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
extdof = find(P(1,:).^2 + P(2,:).^2 > (1 - ...
hmax)^2 & P(1,:).^2 + P(2,:).^2 <= 1^2);

%% (2) State transitions per cell

% Rates is a function that returns a 1-by-nmig cell array, where n is the
% number of migration 'potentials': migration rates are proportional to the
% potential gradients. Inputs to are U = cell number, Q = scalar field, QI
% is Ncells-by-Ncapacity-by-ninternal local cell data array,
% P = dof2position map, T = time.
Rates = @(U,Q,QI,P,t){1*Q(:,1)};   % here, only migration pressure

% Drate is a function that returns a 1-by-nmig cell array that
% specifies how respective migration rate in Rates is scaled. First entry
% scales the first migration rate in Rates (here grad(Q(:,1))), etc.
% Uf is the cell number in the voxel to move from, Ut in voxel to move to.
Dscale = 100;  % scaling of the movement rates (100=>~3 min to move 1 diam.
% when r0 is full of doubles)
Drate = @(Uf,Ut,Q,QI,P,t){Dscale.*(Uf==1).*(Ut<1)+Dscale.*(Uf==2).*(Ut<2)};
                       % Note, solver considers only migration from
                       % bdof_m and sdof_m!

%% (3) Form population

% initial small population
ii1 =  find(P(1,:).^2 + P(2,:).^2 <= R0^2);   % alive cells

% U is Ntype-by-Ncells sparse vector, representing the cell population
U = zeros(ntypes, Nvoxels^2);
U(1,ii1) = 1; % singly occupied

%% (4) "outer" URDME-struct
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);                    % construct D matrix, etc.

r0 = P(1,:).^2 + P(2,:).^2 <= 0.2^2; % circle of proliferation
% 'UL' means cell of phenotype L.
% Define all reaction events first, and after that quantities:
umod = rparse(umod, { ...
              'U1 > (sd==2)*mu_prol*(U1==1)*grow > U1+U1', ...% prolif.
              'Q1 > (U1>1) > Q1+Q1'}, ...
              {'U1' 'Q1'}, ...
              {'mu_prol', mu_prol, 'grow', 'gdata_time'}, ...
              'basic_test_outer');
umod.data_time = [0 1/mu_prol*5];   % ~the time for population to grow X5
umod.gdata_time = [1 0];            % boolean (on/off)
umod.u0 = [U; zeros(1,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(r0) = 2;                                % only proliferation in r0
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps

% Generate internal initial state data.
% scaling and cell/voxel volume
% note: concentrations with an estimated spread
conc = hes1_conc(numel(ii1)); % given in micromolar (uM)
avogadro = 6.022*1e23;
% 1 VOL = 1 um^3 - scaling to get volume in litres
cell_vol = VOL*1e-15;
mumod.vol = repmat(cell_vol,1,numel(umod.vol));
mumod.sd = ones(1,numel(umod.sd));

% fetch the reactions
mumod = hes1umod_dlcm(mumod);

% go from initially given concentrations to # molecules:
%   moles [mol] = molar concentration [M] * volume [L]
%   # molecules = moles [mol] * avogadro's number [#/mol]
scaling_factor = zeros(size(conc,1), Nvoxels^2);
scaling_factor(:,ii1) = conc*1e-6*cell_vol*avogadro;

% initial values (close to given concentrations but in # molecules)
u0 = mumod.u0.*scaling_factor; % scaling_factor([1 1:5],:) before
% "stochastic rounding":
mumod.u0 = floor(u0)+(rand(size(u0)) < u0-floor(u0));

% random (0,Dmax or Nmax) initial delta & notch values in each cell:
nstates = 5;                 % nr of distinct internal states

% celldata is an Ntypes*2-by-Ncells array of cell state data -> mumod.u0.
% Here containing Delta and Notch values for each cell
celldata = zeros(nstates*2, Nvoxels^2);
adof = find(sum(U,1)>0); sdof = find(sum(U,1)>1);
celldata(1:2:end, :) = mumod.u0; % simply extending the initial state
% into including empty 2nd cells as well

%% (5) Internal States using URDME SSA

% Maximum micro-time step (during which ldata is fixed) by "5%-norm rule"
maxdt_fun = @(U, Q, QI, ldata_fun)(0.1); %

% Define reaction rates ldata[0:nldata] (held constant during maxdt)
% Function for ldata input that returns a 1-by-nldata cell array, with
% ldata_fun(...){l} containing ldata[l] for l = 1:nldata.
% Here the function return the avg. neighboring Delta activity
ldata_fun = @(U, Q, QI, P, Ne){sum(Ne*QI(:,:,3),2)./max(Ne*U,1)};

% Compile umod as usual, but need to include diffusion artificially
% in rates, since SSA solver can't include it currently...
% Construct, parse, and compile ssa solver fully outside dlcm
Dexpr = cell(2*nstates,1); Dexpr(1:2*nstates) = {1};
mumod_tmp = pde2urdme(P,T,Dexpr);   % get D
mumod.D = mumod_tmp.D;
mumod.D = sparse(size(mumod.D,1), size(mumod.D,2));
mumod.pde = mumod_tmp.pde;    % vol and sd already defined
mumod.seed = seed;
Nlive = numel(ii1);                     % nr of live cells
Nspec = nstates*2;
mumod.u0 = [];  % dirty solution ('forget' mumod.u0 defined above
mumod.u0(1:Nspec,:) = celldata;         % internal states
mumod.tspan = [0 1];                    % compile dummy
mumod = urdme(mumod, 'solve', 0, 'solver','ssa');       % parse & compile

% finally load essentials and mumod with internal specifics into umod
umod = dlcm2urdme(umod, P, gradquotient, [],[],[], 'Rates', Rates, ...
  'Drate', Drate, 'mumod', mumod, 'ldata_fun', ldata_fun, ...
  'maxdt_fun', maxdt_fun, 'internal_state', internal_state);

%% (6) solve

% Solver automatically neglects internal states (other than types) and
% curvature evaluation if these objects are not defined in solverargs
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', seed);
