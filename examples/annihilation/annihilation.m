% URDME example file. Species A and B start at opposite boundaries (in 1D) 
% and annihilate when they meet. 
%
% This file examplifies how inline propensities are used, and how on can
% set up a 1D simulation. 
%
% B. Drawert, 2012. 


function umod = annihilation(umod)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D_const = 1;
R_const = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mspecies   = 2;
Mreactions = 3;

%connectivity
G=sparse(zeros(Mreactions,Mspecies+Mreactions));
G(3,:)=[1 1 1 1 1];
umod.G = G;

s_A=1;
s_B=2;

% Stoichiometry matrix.
N=zeros(Mspecies,Mreactions);
% 0->A
N(s_A,1)=1;
% 0->B
N(s_B,2)=1;
% A+B->0
N([s_A s_B],3)=[-1 -1];
%%%
umod.N = sparse(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_react=R_const*1.0;
k_creat=R_const*200.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
umod.M1=3;
umod.K = zeros(3,umod.M1);
umod.I = zeros(3,umod.M1);
umod.S = sparse(zeros(3,umod.M1));
% Reaction 1: 0->A (sd == 2) return data[0]*k_creat*vol;
umod.K(3,1) = k_creat; %first index=3 for zero-order reactions
umod.S([1 3],1)=[1;1]; %disable in subdomain 1 and 3
% Reaction 2: 0->B (sd == 3) return data[0]*k_creat*vol;
umod.K(3,2) = k_creat; %first index=3 for zero-order reactions
umod.S([1 2],2)=[1;1]; %disable in subdomain 1 and 2
% Reaction 3: A+B->0 
umod.K(1,3) = k_react; %first index=1 for second-order reactions
umod.I([1 2],3) = [s_A;s_B]; %Bi-molecular reaction reagents 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumParameters=1;
umod.parameters=zeros(Mreactions,NumParameters);
umod.parameters(:,1) = [k_creat;k_creat;k_react];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get # of mesh points
L = 10;

% PB: varargin has to be only an umod object , so you have to include your
% parameters into it
Npoints=20;

x = linspace(0,L,Npoints);
dx = L/Npoints;

umod.mesh.p = zeros(3,Npoints);
umod.mesh.p(3,:) = x;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ncells = numel(x);
Ndofs  = Mspecies*Ncells;


DC = ones(Mspecies,1)*D_const;

% Assemble interior cells
D = spalloc(Ndofs,Ndofs, 3*Ndofs);
h = dx*dx;

for i = 2:Ncells-1
   for j=1:Mspecies
      dof = (i-1)*Mspecies+j;
      D(dof,dof) = -2*DC(j)/h;
      D(dof,dof+Mspecies) = DC(j)/h;
      D(dof,dof-Mspecies) = DC(j)/h;
   end
end
% Reflecting Boundary conditions
for j=1:Mspecies
  dof = j;
  D(dof,dof) = -1*DC(j)/h;
  D(dof,dof+Mspecies) = DC(j)/h;
  dof = (Ncells-1)*Mspecies+j;
  D(dof,dof) = -1*DC(j)/h;
  D(dof,dof-Mspecies) = DC(j)/h;
end

% Compress
umod.D   = sparse(D);
umod.vol = (dx^3)*ones(Ncells,1);
% Three subdomains
umod.sd  = ones(1,Ncells);
umod.sd(1)=2;
umod.sd(end)=3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Initial Conditions
umod.u0 = zeros(Mspecies,Ncells);

% Setup  sample times (every second)
umod.tspan = 0:1:200;

% empty data vec
umod.data = zeros(0,Ncells);