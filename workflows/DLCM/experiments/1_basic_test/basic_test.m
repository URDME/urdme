%BASIC_TEST Basic test of DLCM relaxation.
%   Relaxation to equilibrium: here the cells start distributed within
%   an internal square region filled with 2 cells (red) per voxel. The
%   system then relaxes to equilibrium by, at any given point in time,
%   assuming quasi steady-state and solving an equation for the 'cellular
%   pressure'.

%   Outline of method: at any given point in time we assume quasi
%   steady-state and solve an equation for the 'cellular pressure' in
%   the form of a Laplacian with source terms. Cells can move only
%     -when they have empty voxels next to them (i.e. at boundary
%     points), and here the gradient of the pressure in that direction
%     is understood as a rate per unit of time to change position,
%     -or in general, when a neighbor voxel is less populated than the
%     current one, then the (positive) gradient of the pressure in
%     that direction is again understood as a rate per unit of time to
%     move.
%
%   This stochastic process is simulated in continuous time in the
%   form of a Markov chain.
%
%   The algorithm thus consists of two steps:
%     (1) assume quasi equlibrium (no cells make any large movements,
%     only small movements about each center of mass) - solve the
%     pressure equation with sources at each cell position where there
%     are > 1 cells,
%     (2) all rates determined in this way now imply a cell which can
%     move - find out which ones moves first, and move it.

% E. Blom 2024-12-17 (major revision, using URDME's DLCM solver)
% S. Engblom 2017-12-20 (revision, more cleanup)
% S. Engblom 2017-08-29 (revision, cleanup)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-25 (seriously finalized the physics)
% S. Engblom 2016-12-20 (finalized the physics)
% S. Engblom 2016-12-09 (hexagonal mesh)
% S. Engblom 2016-12-08 (notes on thinning)
% S. Engblom 2016-12-02 (revision)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

%% (1) geometry
Nvoxels = 41;
mesh_type = 1;  % cartesian mesh

% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 30000;
end
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
U(1,:) = fsparse(ii1(:),1,2,[Nvoxels^2 1]); % doubly occupied

%% (4) "outer" URDME-struct
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);                    % construct D matrix, etc.
% 'UL' means cell of phenotype L.
% Define all reaction events first, and after that quantities:
umod = rparse(umod, { ...
              ...%'U1 > 0.0001*(U1==1) > U1+U1', ...% <- adds proliferation
              'Q1 > (U1>1) > Q1+Q1'}, ...
              {'U1' 'Q1'}, ...
              {}, ...
              'basic_test_outer');
umod.u0 = [full(U); zeros(1,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essentials
umod = dlcm2urdme(umod, P, gradquotient, [], [], [], 'Rates', Rates, ...
  'Drate', Drate);

%% (5) solve

% Solver automatically neglects internal states (other than types) and
% curvature evaluation if these objects are not defined in solverargs
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
