%AVASCULAR_TUMOR Avascular tumor model
%   Small intial population of cells that proliferate based on a
%   nutrient diffusing in from a 'far-away' external source (blood
%   vessels). As the tumor grows, nutrients become scarce, resulting
%   in degradation of cells. Surface tension on the population
%   boundary can keep the tumor boundary circular when nutrients are
%   scarce.

%   Avascular tumour growth: An initial circular population cells (one
%   per voxel) lie in a domain rich in oxygen. Cells consume oxygen at
%   a constant rate, lambda. Cells occupying a voxel with oxygen above
%   kappa_prol can proliferate at a rate mu_prol. Cells occupying
%   voxels with an oxygen concentration below kappa_die can die at a
%   rate mu_die.  Dead cells are represented as a different cell type and
%   can degrade to stop occupying space at a rate mu_deg.
%
%   No oxygen consumption by dying cells.
%   Surface tension on every contour that is 'large enough'
%
%   Permeability: Drate describes the scaling of the migration rate of
%   tumor cells invading occupied or empty voxels.

% E. Blom 2025-04-17 (major revision, using URDME's DLCM solver)
% E. Blom 2023-03-20 (revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

%% (1) geometry and surface tension
Nvoxels = 51;
mesh_type = 2;  % hexagonal mesh

% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 10000;
end
Tres = 100;
% surface tension coefficients
ntypes = 2; % number of cell types: 1: degrading 2: living
% sigma counts type 0 as void:
sigma = zeros(ntypes+1, ntypes+1);  % sigma_ij: type i-1 <-> type j-1
sigma(2:ntypes+1,1) = [0e-4 0e-4];  % zero sigma => tumor breaks
sigma(3:ntypes+1,2) = 0e-4;
sigma = sigma+sigma';               % should be symmetric, zero diagonal

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.

% get boundary dofs, extdof
hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
extdof = find(P(1,:).^2 + P(2,:).^2 > (0.866 - ...
hmax)^2 & P(1,:).^2 + P(2,:).^2 <= 0.866^2);

%% (2) Migration rates per cell

% Rates is a function that returns a 1-by-nmig cell array, where n is the
% number of migration 'potentials': migration rates are proportional to the
% potential gradients. Inputs to are U = cell number, Q = scalar field, QI
% is Ncells-by-Ncapacity-by-ninternal local cell data array,
% P = dof2position map, T = time.
Rates = @(U,Q,QI,P,t){Q(:,1)};

% Drate is a function that returns a 1-by-nmig cell array that
% specifies how respective migration rate in Rates is scaled. First entry
% scales the first migration rate in Rates (here grad(Q(:,1))), etc.
% Uf is the cell number in the voxel to move from, Ut in voxel to move to.
Drate = @(Uf,Ut,Q,QI,P,t){1.*(Uf==1).*(Ut<=1)+...
                         (Uf==2).*(25.*(Ut==1)+1.*(Ut==0))};

%% (3) Form population

% initial small population
ii1 = find((P(1,:)).^2 + (P(2,:)).^2 <= 0.05^2);   % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(2,:) = fsparse(ii1(:),1,1,[Nvoxels^2 1]); % no degrading cells initially

%% (4) "outer" URDME-struct
% Living cells can become degrading cells that vanish over time, and living
% cells can proliferate. Living cell events depend on oxygen levels.

% Phenotype switching rates
mu_die = 0.1;             % rates of birth, death, and degradation
mu_deg = mu_die*0.1;
mu_prol = 0.2;
kappa_prol = 0.955;      % oxygen thresholds for birth & death
kappa_die = 0.945;

nquants = 2; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% 'UL' means cell of phenotype L.
% Define all reaction events first, after that quantities:
umod = rparse(umod, ...
              {'U1 > mu_deg*(U1>0) > @', ...        % degrade, birth, death
              'U2 > mu_prol*(Q2>kappa_prol)*(U2==1) > U2+U2', ...
              'U2 > mu_die*(Q2<kappa_die)*(U2>0) > U1', ...
              'Q1 > (U2>1) > Q1+Q1', ...            % pressure (y-l bc)
              'Q2 > sd > 0 ? -U2 : 1 > Q2+Q2'}, ... % set BC at sd == 0
              {'U1' 'U2' 'Q1' 'Q2'}, ...
              {'mu_deg' mu_deg, ...
              'mu_prol' mu_prol, 'mu_die' mu_die, ...
              'kappa_prol', kappa_prol, 'kappa_die', kappa_die}, ...
              'avascular_tumor_outer');
umod.u0 = [full(U); zeros(2,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                              % sd encodes bndry dofs
umod.tspan = linspace(0,Tend,Tres);               % simulation time steps
% load essential objects + curvature operators
umod = dlcm2urdme(umod, P, gradquotient, T, sigma, [], 'Rates', Rates, ...
  'Drate', Drate);

%% (5) solve
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
