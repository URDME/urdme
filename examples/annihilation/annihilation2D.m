function [umod,G,P,E,T] = annihilation2D
% Annihilation 2D example.

% S. Engblom 2017-02-21

% scaling of reactions and diffusions
R_const = 1;
D_const = 0.1;

% build the PDE Toolbox model

% (1) create the geometry, here composed of three circles
C1 = [1 0 0 1]';
C2 = [1 -0.4 0 0.2]';
C3 = [1 0.4 0 0.2]';
gd = [C1 C2 C3];
sf = 'C1+C2+C3';
ns = char('C1','C2','C3')';
G = decsg(gd,sf,ns);

% (2) create the mesh
[P,E,T] = initmesh(G,'hmax',0.075);

% (3) assemble the diffusion part
umod = pde2urdme(P,T,{D_const D_const});

% this "rounding" will include only the nodes truly inside C2, C3:
umod.sd(umod.sd < 1.5) = 1;
umod.sd(umod.sd > 2.5) = 3;
umod.sd = round(umod.sd);

% reaction topology:
%
% @ --> A (sd == 1), @ --> B (sd == 3), A+B --> @ (everywhere)
umod.N = sparse([1 0 -1; ...
                 0 1 -1]);
umod.G = sparse([0 0 0 0 0; ...
                 0 0 0 0 0; ...
                 1 1 1 1 1]);

% inline propensities
k_react = R_const*1.0;
k_creat = R_const*200.0;

K = [0 0 k_react; ...
     0 0 0; ...
     k_creat k_creat 0];
I = [1 1 1; ...
     1 1 2; ...
     1 1 1];
S = [2 1 0; ...
     3 2 0];
S = sparse(S);

umod.solverargs = {'K' K 'I' I 'S' S};

% initial conditions
Mspecies = size(umod.N,1);
Ncells = numel(umod.vol);
umod.u0 = zeros(Mspecies,Ncells);

% sample times
umod.tspan = 0:100-1;

% empty data
umod.ldata = zeros(0,Ncells);
umod.gdata = [];
