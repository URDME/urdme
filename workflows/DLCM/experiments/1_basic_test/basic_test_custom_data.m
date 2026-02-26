%BASIC_TEST_CUSTOM_DATA Basic test of DLCM relaxation with custom data.
%   Similar to basic_test, but tests custom ldata and gdata
%   functionality using (dummy) internal states.

% E. Blom 2026-01-07

rng(123) % for birth rates in ldata

%% (1) geometry
Nvoxels = 41;
mesh_type = 1;  % cartesian mesh

% Simulate to Tend and save states at Tres intervals
Tend = 100;
Tres = 100;
ntypes = 1; % number of cell types: living cells

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.

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
Rates = @(U,Q,QI,P,t){Q(:,1)};   % here, only migration pressure

% Drate is a function that returns a 1-by-nmig cell array that
% specifies how respective migration rate in Rates is scaled. First entry
% scales the first migration rate in Rates (here grad(Q(:,1))), etc.
% Uf is the cell number in the voxel to move from, Ut in voxel to move to.
Drate = @(Uf,Ut,Q,QI,P,t){1.*(Uf==1).*(Ut==0)+1.*(Uf==2).*(Ut<2)};
                       % Note, solver considers only migration from
                       % bdof_m and sdof_m!

%% (3) Form population

% initial small population
ii1 =  find(abs(P(1,:)) < 0.4 & abs(P(2,:)) <= 0.4);   % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U = zeros(ntypes, Nvoxels^2);
U(1,ii1) = 2; % doubly occupied

%% (4) "outer" URDME-struct
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);                    % construct D matrix, etc.
% 'UL' means cell of phenotype L.
% Define all reaction events first, and after that quantities:
umod = rparse(umod, { ...
              'U1 > mu_prol*(U1==1) > U1+U1', ...
              'U1 > ldata[4]*(U1>0) > @', ... % access the data like this...
              'Q1 > p_source*(U1>1) > Q1+Q1'}, ...
              {'U1' 'Q1'}, ...
              {'mu_prol', 'gdata', 'p_source', 'gdata'}, ... % or like this
              'basic_test_outer');
umod.u0 = [full(U); zeros(1,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
% custom ldata (2nd row is dummy):
% Note that ldata [4] corresponds to first array below ([0] & [1] is the
% cell's position map for x and y and [2] & [3] corresponds to the internal
% state (one per possible cell in the voxel), which is set to zero below.
umod.ldata = [1*0.001*rand(1,Nvoxels^2); exp(1)*ones(1,Nvoxels^2)];
% custom gdata:
umod.gdata = [0.0002 1];
umod.tspan = linspace(0,Tend,Tres);             % time steps

% add dummy internal states to check that solver handles them correctly
% when also using custom ldata
nstates = 1;
% required to pass internal states through mumod:
mumod.seed = 1;
mumod.u0(1:2*nstates,:) = 0;         % mumod holds internal states

% load essentials
umod = dlcm2urdme(umod, P, gradquotient, [], [], [], 'Rates', Rates, ...
  'Drate', Drate, 'mumod', mumod);

%% (5) solve

% Solver automatically neglects internal states (other than types) and
% curvature evaluation if these objects are not defined in solverargs
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
