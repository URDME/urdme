%BASIC_TEST_3D Basic test of DLCM relaxation in 3D.
%   3D variant of BASIC_TEST. Initially cubic population with singly
%   occupied voxels and a sphere with doubly occupied voxels half-inside it.
%   Uses 3D curvature operators for population surface tension.

% E. Blom 2024-12-17 (major revision, using DLCM's URDME solver)
% S. Engblom 2017-08-29 (revision, cleanup)
% S. Engblom 2016-11-13 (revision, 3D version)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

%% (1) geometry
% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 2000000;
end
Tres = 100;
ntypes = 1; % number of cell types: living cells
sigma = zeros(ntypes+1, ntypes+1);    % sigma_ij: type i-1 <-> type j-1
sigma(2:ntypes+1,1) = [1e-5];         % 1 with void = 0
sigma = sigma+sigma';

% construct 3D sphere
hmax = 0.08;
gm = fegeometry(multisphere(1));
gm = generateMesh(gm, 'GeometricOrder', 'linear', Hmax=hmax);
[P,E,T] = meshToPet(gm.Mesh);
Nvoxels = size(P,2);

% calculate individual tetrahedron edge areas
%t = 220;     % arbitrary edge index
%v = [P(1, T(1:3,t)); P(2, T(1:3,t)); P(3, T(1:3,t)) ];
%0.5*sqrt(sum(v(:,1).*v(:,1))*sum(v(:,2).*v(:,2)) - sum(v(:,1).*v(:,2))^2)
edge_area = 0.02;    % roughly 0.02 evaluated manually using above...
gradquotient = edge_area*2*(gm.Mesh.MaxElementSize+gm.Mesh.MinElementSize);

% define boundary nodes ('extdofs' in the framework)
% crude approach to finding the outermost nodes on the sphere
hmin = gm.Mesh.MinElementSize;
extdof = find(P(1,:).^2 + P(2,:).^2 + P(3,:).^2 > (1 - ...
  hmin)^2 & P(1,:).^2 + P(2,:).^2 + P(3,:).^2 <= 1);

%% (2) State transitions per cell

% Rates is a function that returns a
% 1-by-nevents+1 cell array: the rate of event in N(:,i) is defined by
% Rates(.){i}. The final cell Rates(.){end} defines the migration rate
% functional, whose gradient gives the migration rate. Inputs to qrhs(.)
% are U = cell number, Q = scalar field, QI is Ncells-by-Ncapacity-by-
% -ninternal local cell data array, P = dof2position map, T = time.
Rates = @(U,Q,QI,P,t){Q(:,1)};   % here, only migration pressure

% Drate is a function that returns a 1-by-Ncapacity cell array that
% specifies how migration rates are scaled. Uf is the cell number in the
% voxel to move from, Ut in voxel to move to.
Drate = @(Uf,Ut,Q,QI,P,t){1.*(Uf==1).*(Ut<=1)+1.*(Uf==2).*(Ut<2)};
                       % Note, solver considers only migration from
                       % bdof_m and sdof_m!

%% (3) Form population

% initial small population
ii1 =  find(abs(P(1,:)) < 0.35 & abs(P(2,:)) <= 0.35 & abs(P(3,:)) <= 0.35);
ii1_b =  find(sqrt((P(1,:)-0.4).^2 + (P(2,:)).^2 + (P(3,:)).^2) <= 0.3);
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(1,:) = fsparse(ii1(:),1,1,[Nvoxels 1]); % singly occupied cube
U(1,ii1_b) = 2;                           % doubly occupied sphere closeby

%% (4) "outer" URDME-struct
nquants = 1; % number of field states pressure and nutrient
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);                    % construct D matrix, etc.
% 'UL' means cell of phenotype L.
% !Define all reaction events first, after that quantities!
umod = rparse(umod, ...
              {'Q1 > (U1>1) > Q1+Q1'}, ...
              {'U1' 'Q1'}, ...
              {}, ...
              'basic_test_3D_outer');
umod.u0 = [full(U); zeros(1,Nvoxels)];
umod.sd = ones(1,Nvoxels);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essentials
umod = dlcm2urdme(umod, P, gradquotient, T, sigma, [], 'Rates', Rates, ...
  'Drate', Drate);

%% (6) solve

% Solver automatically neglects internal states (other than types) and
% curvature evaluation if these objects are not defined in solverargs
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
