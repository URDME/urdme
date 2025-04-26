%CHEMOTAXIS_MODLA Chemotaxis model using custom Laplacian.

% E. Blom 2024-12-02

%% (1) geometry
% Simulate to Tend and save states at Tres intervals
if ~exist('Tend', 'var')
  Tend = 2000000;
end
Tres = 100;
ntypes = 1; % number of cell types: living cells

Nvoxels = 101;
mesh_type = 2;  % hexagonal mesh

% Laplacian and neighbor matrix assembled manually
[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh - for visualisation

% (extdofs not required here, the only BCs are on pupulation boundary)

%% (2) assemble custom Laplacian and neighbor matrix:

% gradient of the SLIT chemical
% parameters for the diffusion and degradation of the chemical
Q = 2; k = 0.1; Ds = 50; lambda = sqrt(k/Ds);
sx = Q/(2*Ds*lambda).*exp(-lambda.*(1+P(1,:)));

[L,dM,~,~,~] = dt_operators(P,T);

% manually assemble convection matrix (will use its transpose)
chi = -2e4;           % diffusion factor (large factor needed)
bx = -chi*lambda*sx;  % derivative in x-direction
by = 0*sx;            % derivative in y-direction (SLIT is constant in y
np=size(P,2);
nt=size(T,2);
C=sparse(np,np);
for i=1:nt
  % Assembly from M.G Larson, F. Bengzon, "The Finite Element Method:
  % Theory, Implementation, and Applications" (2013)
  loc2glb=T(1:3,i);
  x=P(1,loc2glb);
  y=P(2,loc2glb);
  [area,b,c]=HatGradients(x,y);
  bxmid=mean(bx(loc2glb));
  bymid=mean(by(loc2glb));
  CK=ones(3,1)*(bxmid*b+bymid*c)'*area/3;
  C(loc2glb,loc2glb)=C(loc2glb,loc2glb)+CK;
end

% advection-diffusion operator
C = sparse(diag(1./dM)*C');
A = (L + C);
% neighbor matrix
N = (A-diag(diag(A))~=0);

%% (3) State transitions per cell

% Here we assume a variant of chemotaxis where the pressure and signalling
% substance are coupled as one quantity determining the migration
Rates = @(U,Q,QI,P,t){Q(:,1)}; % chemical-pressure coupling

% Cells move up the chemical gradient by some innate sensitivity
Drate = @(Uf,Ut,Q,QI,P,t){((Uf==1) + 1*(Uf==2)).*(Ut==0)};

%% (4) Form population

% initial small population
ii1 = find((P(1,:)).^2 + (P(2,:)).^2 <= 0.25^2); % alive cells
% U is Ntype-by-Ncells sparse vector, representing the cell population
U(1,:) = fsparse(ii1(:),1,2,[Nvoxels^2 1]);

%% (5) URDME struct
nquants = 1; % number of field states, only pressure here
Dexpr = cell(1,nquants+ntypes);
Dexpr(:) = {1};
umod = pde2urdme(P,T,Dexpr);
% 'UL' means cell of phenotype L.
% Define all reaction events first, after that quantities
umod = rparse(umod, ...
              {'Q1 > (U1>1) > Q1+Q1'} , ...    % pressure (y-l bc)
              {'U1' 'Q1'}, ...
              {}, ...
              'chemotaxis3D_outer');
umod.u0 = [full(U); zeros(nquants,Nvoxels^2);];
umod.sd = ones(1,Nvoxels^2);
umod.tspan = linspace(0,Tend,Tres);             % time steps
% load essential objects
umod = dlcm2urdme(umod, P, gradquotient, [],[],[], 'Rates', Rates, ...
  'Drate', Drate, 'D', A, 'Ne', N);

%% (6) solve (with custom D and N)
umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123);

return;

% uncomment to do multiple sample runs:
% U_samples = 0;
% Nsamples = 100;
% chemotaxis_modLa
% for i = 1:Nsamples
%   U_samples = U_samples + reshape(sum(umod.U(1:ntypes,:,:),1), ...
%   size(umod.U,2),Tres)./Nsamples;
%   % (last solve is not counted)
%   umod = urdme(umod,'solver','dlcm', 'solve', 1, 'seed', 123+10*i);
%   i
% end

function [area,b,c] = HatGradients(x,y)
%   Basic function to evaluate the constant gradients of the hat
%   functions for Finite Element Methods.
%   From M.G Larson, F. Bengzon, "The Finite Element Method:
%   Theory, Implementation, and Applications" (2013).

area=polyarea(x,y);
b=[y(2)-y(3); y(3)-y(1); y(1)-y(2)]/2/area;
c=[x(3)-x(2); x(1)-x(3); x(2)-x(1)]/2/area;
end
