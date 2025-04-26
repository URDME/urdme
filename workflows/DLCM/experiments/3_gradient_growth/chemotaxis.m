%CHEMOTAXIS Chemotaxis models.
%   Model_type 1: Darcy's law + chemotaxis towards attracting chemical
%   Model_type 2: Diffusion +   --||--
%   Model_type 3: Darcy's law + chemo-repellant that cells consume

% E. Blom 2024-12-02 (major revision, URDME solver & additional models)
% S. Engblom 2017-12-28 (minor revision)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

%% (1) geometry
% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 30000;
end
Tres = 100;
ntypes = 1;     % number of cell types: living cells

Nvoxels = 101;
mesh_type = 2;  % hexagonal mesh

if ~exist('model_type', 'var')
  model_type = 1; % [1, 2, 3]
end

[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - only needed for vis.

% get boundary dofs, extdof (Neumann on boundaries not in extdof!)
hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
keep = find(P(1,E(1,:)) < -0.99 | P(1,E(1,:)) > 0.99); % only at x=1 & -1
extdof = unique(E(1,keep));

%% (2) State transitions per cell

% Here we assume a variant of chemotaxis where the pressure and signalling
% substance are decoupled and we take their gradients separately when
% evaluating the cell migration
switch model_type
  case 1
    % attractant
    Rates = @(U,Q,QI,P,t){Q(:,1), -Q(:,2)};
  case 2
    % (diffusion instead of Darcy's law:)
    Rates = @(U,Q,QI,P,t){0.001*U, -Q(:,2)};
  case 3
    % repellant
    Rates = @(U,Q,QI,P,t){Q(:,1), Q(:,2)};
end

% cells move up the chemical gradient by some innate sensitivity QI(:,1,3),
% here constant but is easily adjusted for more detail.
Drate = @(Uf,Ut,Q,QI,P,t){(Uf==1).*(Ut==0)+(Uf==2).*(Ut<=1), ...
                      QI(:,1,3).*((Uf==1)+(Uf==2)).*(Ut<=1)};

%% (3) Form population

% initial small population
ii1 = find((P(1,:)).^2 + (P(2,:)).^2 <= 0.25^2); % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(1,:) = fsparse(ii1(:),1,1,[Nvoxels^2 1]);

% Define internal states, here chemotactic sensitivity value = 1
% celldata is an Ncells-by-Ncapacity-by-ninternal array with internal state
% data corresponding to the format of QI in InternalRates defined above
celldata = zeros(ntypes*2, Nvoxels^2);
adof = find(U(1,:)>0); sdof = find(U(1,:)>1);
celldata(1, adof) = ones(numel(adof), 1);   % first cell in voxel
celldata(2, sdof) = ones(numel(sdof), 1);   % second cell in voxel

nstates = 1;
% required to pass internal states through mumod:
mumod.seed = 1;
Nspec = 2*nstates;
mumod.u0(1:Nspec,:) = celldata;         % internal states

%% (4) URDME struct
nquants = 2; % number of field states pressure and chemical signal
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% model 3: signalling consumption and attracting chemical (to the left)
r2 = 'Q2 > sd>0 ? -0.1*U1: 0.1*(ldata[0]>0.978) > Q2+Q2';
if model_type <= 2
  % repelling chemical (from the right)
  r2 = 'Q2 > sd>0 ? 0 : 0.1*(ldata[0]<-0.978) > Q2+Q2';
end
% 'UL' means cell of phenotype L.
% !Define all reaction events first, after that quantities!
umod = rparse(umod, ...
              {'Q1 > (U1>1) > Q1+Q1', ... % pressure (y-l bc)
              r2}, ...
              {'U1' 'Q1' 'Q2'}, ...
              {}, ...
              'chemotaxis_outer');
umod.u0 = [full(U); zeros(2,Nvoxels^2);];
umod.sd = ones(1,Nvoxels^2);
umod.sd(extdof) = 0;                            % sd encodes boundary dofs
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essentials and internal dynamics in mumod
umod = dlcm2urdme(umod, P, gradquotient, [],[],[], 'Rates', Rates, ...
  'Drate', Drate, 'mumod', mumod);

%% (6) solve
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);

return;

%% do multiple sample runs
U_samples = 0;
Nsamples = 100;
chemotaxis
for i = 1:Nsamples
  U_samples = U_samples + reshape(sum(umod.U(1:ntypes,:,:),1), ...
  size(umod.U,2),Tres)./Nsamples;
  % last solve is not counted...
  umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123+10*i);
  i
end
