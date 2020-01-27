% Annihilation 2D example.
%   Species A and B are created in different subdomains and annihilate
%   when they meet. This example shows how one can set up a simple 2D
%   simulation using the URDME interface to PDE Toolbox.

% S. Engblom 2019-11-27 (Revision, model augmentation using rparse_inline)
% S. Engblom 2017-05-09 (Revision, rparse_inline)
% S. Engblom 2017-02-21

%% (1) geometry and diffusion operator

% scaling of reactions and diffusions
R_const = 1;
D_const = 0.1;

% build the PDE Toolbox model

% create the geometry, here composed of three circles
C1 = [1 0 0 1]';
C2 = [1 -0.4 0 0.2]';
C3 = [1 0.4 0 0.2]';
gd = [C1 C2 C3];
sf = 'C1+C2+C3';
ns = char('C1','C2','C3')';
G = decsg(gd,sf,ns);

% create the mesh
[P,E,T] = initmesh(G,'hmax',0.075);

% assemble the diffusion part
umod = pde2urdme(P,T,{D_const D_const});

% this "rounding" will include only the nodes truly inside C2, C3:
umod.sd(umod.sd < 1.5) = 1;
umod.sd(umod.sd > 2.5) = 3;
umod.sd = round(umod.sd);

%% (2) reactions

% @ --> A (sd == 1), @ --> B (sd == 3), A+B --> @ (everywhere)
k_creat = R_const*200.0;
k_react = R_const*1.0;
umod = rparse_inline(umod, ...
    {'@ > k_creat > A', ...
     '@ > k_creat > B', ...
     'A+B > k_react > @'}, ...
    {'A' 'B'},{'k_creat' k_creat 'k_react' k_react});

% specify in what subdomains reactions are turned off
S = [2 1 0; ...
     3 2 0];
S = sparse(S);
umod.inline_propensities.S = S;

%% (3) run the model

% initial conditions
Mspecies = size(umod.N,1);
Ncells = numel(umod.vol);
umod.u0 = zeros(Mspecies,Ncells);

% simulate
umod = urdme(umod,'tspan',0:1000-1,'seed',20170221);
umod = urdme2pde(umod);

% compute temporal means
mean_A = mean(umod.pde.U(1,:,end/2:end),3);
mean_B = mean(umod.pde.U(2,:,end/2:end),3);

% visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf,
  pdegplot(G,'subdomainlabels','on'), axis equal

  figure(2), clf,
  pdemesh(P,E,T), axis equal
  
  figure(3), clf
  pdesurf(umod.pde.P,umod.pde.T,mean_A');
  title('Mean A');
  
  figure(4), clf
  pdesurf(umod.pde.P,umod.pde.T,mean_B');
  title('Mean B')
end
