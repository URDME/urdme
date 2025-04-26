%CHEMOTAXIS3D Chemotaxis model in three spatial dimensions.
%   Cells emit a chemical signal that attracts other cells. A fraction of
%   the cells have a higher sensitivity to the signal, thus aggregating.

% E. Blom 2024-12-02

rng(123) % for initial states of chemotactic sensitivity

%% (1) geometry
% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 100000;
end
Tres = 100;
ntypes = 1; % number of cell types: living cells

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

% define boundary nodes ('extdofs' in the framework) by
% finding the outermost nodes on the sphere:
hmin = gm.Mesh.MinElementSize;
extdof = find(P(1,:).^2 + P(2,:).^2 + P(3,:).^2 > (1 - ...
  hmin)^2 & P(1,:).^2 + P(2,:).^2 + P(3,:).^2 <= 1);

% uncomment to try the 2D version:
% Nvoxels = 101;
% mesh_type = 2;  % cartesian mesh
%
% [P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
% [V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.
%
% % get boundary dofs, extdof
% hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
% extdof = find(P(1,:).^2 + P(2,:).^2 > (0.866 - ...
% hmax)^2 & P(1,:).^2 + P(2,:).^2 <= 0.866^2);
% Nvoxels = Nvoxels^2;

%% (2) State transitions per cell

% Here we assume a variant of chemotaxis where the pressure and signalling
% substance are decoupled and we take their gradients separately when
% evaluating the cell migration
Rates = @(U,Q,QI,P,t){Q(:,1), ...   % pressure
    -Q(:,2)};                       % chemical signal

% i) Pressure migration scaling,
% ii) Cells move up the chemical gradient by their innate sensitivity:
Drate = @(Uf,Ut,Q,QI,P,t){1.*((Uf==1)+(Uf==2)).*(Ut<=1); ...  % i)
                      QI(:,1,3).*((Uf==1)+(Uf==2)).*(Ut<=1)}; % ii)

%% (3) Form population

% initial small population
ii1 = find((P(1,:)).^2 + (P(2,:)).^2 + (P(3,:)).^2 <= 0.4^2); % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(1,:) = fsparse(ii1(:),1,1,[Nvoxels 1]);

% Define internal states, here random chemotactic sensitivity values.
% Celldata is a 2*ninternal-by-Ncells array with the initial internal states
nstates = 1;                                % number of internal states
celldata = zeros(nstates*2, Nvoxels);
adof = find(U(1,:)>0); sdof = find(U(1,:)>1);
celldata(1, adof) = rand(numel(adof), 1);   % first cell in voxel
celldata(2, sdof) = rand(numel(sdof), 1);   % second cell in voxel
celldata(1,celldata(1,:)>0.75) = 5;         % extra sensitive cells
celldata(2,celldata(2,:)>0.75) = 5;         % extra sensitive cells

% required to pass internal states through mumod:
mumod.seed = 1;
mumod.u0(1:2*nstates,:) = celldata;         % mumod holds internal states

%% (4) URDME struct
nquants = 2; % number of field states pressure and chemical signal
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% 'UL' means cell of phenotype L.
% Define all reaction events first, after that quantities
umod = rparse(umod, ...
              {'Q1 > 10*(U1>1) > Q1+Q1', ...      % pressure
              'Q2 > U1-Q2 > Q2+Q2'}, ...          % chemical signal
              {'U1' 'Q1' 'Q2'}, ...
              {}, ...
              'chemotaxis3D_outer');
umod.u0 = [full(U); zeros(2,Nvoxels);];
umod.sd = ones(1,Nvoxels);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essential objects
umod = dlcm2urdme(umod, P, gradquotient, [], [], [], 'Rates', Rates, ...
  'Drate', Drate, 'mumod', mumod);
%% (6) solve
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);
