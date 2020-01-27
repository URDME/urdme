% Annihilation example.
%   Species A and B start at opposite boundaries (in 1D) and
%   annihilate when they meet. This example shows how inline
%   propensities are used, and how one can set up a simple 1D
%   simulation.

% S. Engblom 2019-11-27 (Revision, model creation using rparse_inline)
% S. Engblom 2017-05-09 (Revision, rparse_inline)
% S. Engblom 2017-02-18 (Major revision, URDME 1.3, Comsol 5)
% B. Drawert 2012

%% (1) reactions

% scaling of reactions and diffusions
R_const = 1;
D_const = 1;

% @ --> A (sd == 2), @ --> B (sd == 3), A+B --> @ (everywhere)
clear umod
k_creat = R_const*200.0;
k_react = R_const*1.0;
umod = rparse_inline([], ...
    {'@ > k_creat > A', ...
     '@ > k_creat > B', ...
     'A+B > k_react > @'}, ...
    {'A' 'B'},{'k_creat' k_creat 'k_react' k_react});

% specify in what subdomains reactions are turned off
S = [1 1 0; ...
     3 2 0];
S = sparse(S);
umod.inline_propensities.S = S; 

%% (2) diffusion

% length of domain
L = 10;

% discretization
Ncells = 20;
x = linspace(0,L,Ncells);
dx = x(2)-x(1);
umod.private.mesh = x; % (later used to plot)
umod.vol = dx^3*ones(Ncells,1);

% 1D diffusion operator
e = D_const*ones(Ncells,1)/dx^2;
D = spdiags([e -2*e e],-1:1,Ncells,Ncells);

% reflecting boundary
D(1,1) = -e(1);
D(end,end) = -e(end);

% URDME ordering of species
umod.D = kron(D,eye(2));

% three subdomains
umod.sd = ones(1,Ncells);
umod.sd([1 end]) = [2 3];

%% (3) run the model

% initial conditions
Mspecies = size(umod.N,1);
umod.u0 = zeros(Mspecies,Ncells);

% simulate
umod = urdme(umod,'tspan',0:1000-1,'seed',20170219);

%% (4) postprocessing

% compute means
mean_A = mean(umod.U(1:2:end,end/2:end),2)./umod.vol;
mean_B = mean(umod.U(2:2:end,end/2:end),2)./umod.vol;

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf
  plot(umod.private.mesh,mean_A,'-sb',umod.private.mesh,mean_B,'-dr');
  legend('A','B');
  xlabel('x');
end
