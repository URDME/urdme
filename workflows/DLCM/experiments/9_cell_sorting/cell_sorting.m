%CELL_SORTING Cell sorting model.
%   Two different cell types in initially random distribution by
%   Young-Laplace pressure drops between interfaces. ~0.5-1hr runtime

% E. Blom 2024-11-19

%% (1) geometry and surface tension
Nvoxels = 71;
mesh_type = 2;  % hexagonal mesh
% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 200000;
end
Tres = 100;
% surface tension coefficients
ntypes = 3; % number of cell types
sigma = zeros(ntypes+1, ntypes+1);      % sigma_ij: type i-1 <-> type j-1
sigma(2:ntypes+1,1) = [2e-4 4e-4 6e-4];      % 1,2,3  with void = 0
sigma(3:ntypes+1,2) = [1e-4 1e-4];           % 2,3    with type 1
sigma(4:ntypes+1,3) = [1e-4];                % 3      with type 2
sigma = sigma+sigma';                   % should be symmetric, zero diagonal

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.

% get boundary dofs, extdof
hmax = 2/(Nvoxels-1);
extdof = find(P(1,:).^2 + P(2,:).^2 > (0.866 - ...
hmax)^2 & P(1,:).^2 + P(2,:).^2 <= 0.866^2);

%% (2) State transitions per cell

% Rates is a function that returns a 1-by-nmig cell array, where n is the
% number of migration 'potentials': migration rates are proportional to the
% potential gradients. Inputs to are U = cell number, Q = scalar field, QI
% is Ncells-by-Ncapacity-by-ninternal local cell data array,
% P = dof2position map, T = time.
Rates = @(U,Q,QI,P,t){Q(:,1)};

% Drate is a function that returns a 1-by-Ncapacity cell array that
% specifies how migration rates are scaled. Uf is the cell number in the
% voxel to move from, Ut in voxel to move to.
Drate = @(Uf,Ut,Q,QI,P,t){1.*(Ut<=1)};

%% (3) Form population

% initial small population
ii1 = find((P(1,:)).^2 + (P(2,:)).^2 <= 0.33^2); % cell type 1
ii2 = ii1;   % cell type 2
ii3 = ii1;   % cell type 3
% scramble initial cells
rng(123);
ii2 = randsample(ii2, round(numel(ii2)*0.33));
ii3 = randsample(ii3, round(numel(ii3)*0.33));
ii1 = setdiff(ii1, ii2);    % non-overlapping types
ii1 = setdiff(ii1, ii3);
ii2 = setdiff(ii2, ii3);
ii = [ii1 ii2 ii3];                               % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(1,:) = fsparse(ii1(:),1,1,[Nvoxels^2 1]);
U(2,:) = fsparse(ii2(:),1,1,[Nvoxels^2 1]);
U(3,:) = fsparse(ii3(:),1,1,[Nvoxels^2 1]);

%% (4) outer URDME struct
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% 'UL' means cell of phenotype L.
% Define all reaction events first, and after that quantities:
warning('off', 'rparse:ghost_species')  % U1-3 necessary, but no reactions
umod = rparse(umod, ...
              {'Q1 > (U1 > 1) > Q1+Q1'}, ...
              {'U1' 'U2' 'U3' 'Q1'}, ...
              {}, ...
              'cell_sorting_outer');
umod.u0 = [full(U); zeros(1,Nvoxels^2)];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essential objects +  curvature operators
umod = dlcm2urdme(umod, P, gradquotient, T, sigma, hmax, 'Rates', Rates, ...
  'Drate', Drate);
warning('on', 'rparse:ghost_species')
%% (5) solve
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
