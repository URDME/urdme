function umod = annihilation
% Annihilation example.

% S. Engblom 2017-02-18 (Major revision, URDME 1.3, Comsol 5)
% B. Drawert 2012

% scaling of reactions and diffusions
R_const = 1;
D_const = 1;

% reaction topology:
%
% @ --> A (sd == 2), @ --> B (sd == 3), A+B --> @ (everywhere)
N = [1 0 -1; ...
     0 1 -1];
umod.N = sparse(N);
G = [0 0 0 0 0; ...
     0 0 0 0 0; ...
     1 1 1 1 1];
umod.G = sparse(G);

% inline propensities
k_react = R_const*1.0;
k_creat = R_const*200.0;

K = [0 0 k_react; ...
     0 0 0; ...
     k_creat k_creat 0];
I = [1 1 1; ...
     1 1 2; ...
     1 1 1];
S = [1 1 0; ...
     3 2 0];
S = sparse(S);

umod.solverargs = {'K' K 'I' I 'S' S};

% length of domain
L = 10;

% discretization
Ncells = 20;
x = linspace(0,L,Ncells);
dx = x(2)-x(1);
umod.private.mesh = x; % (used to plot)
umod.vol = dx^3*ones(Ncells,1);

% 1D diffusion operator
e = D_const*ones(Ncells,1)/dx^2;
D = spdiags([e -2*e e],-1:1,Ncells,Ncells);

% reflecting boundary
D(1,1) = -e(1);
D(end,end) = -e(end);

% URDME ordering of species
rep = reshape([1:Ncells; 1:Ncells],1,[]);
umod.D = D(rep,rep);

% three subdomains
umod.sd = ones(1,Ncells);
umod.sd([1 end]) = [2 3];

% initial conditions
Mspecies = size(umod.N,1);
umod.u0 = zeros(Mspecies,Ncells);

% sample times
umod.tspan = 0:1000-1;

% empty data
umod.ldata = zeros(0,Ncells);
umod.gdata = [];
