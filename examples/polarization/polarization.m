function umod = polarization(umod)
% Model file for the yeast polarization example. 
% 

Mspecies   = 6;
Mreactions = 8;

G=sparse(ones(Mreactions,Mreactions+Mspecies));

umod.G = G;

% Stoichiometry matrix.
N=zeros(Mspecies,Mreactions);
% Species
% R RL G Ga Gbg Gd
% 1 2  3 4  5   6
% R1: 0 -> R
N(1,1)=1;
% R2: R -> 0
N(1,2)=-1;
% R3: R -> RL
N([1 2],3)=[-1,1];
% R4: RL -> R 
N([2 1],4)=[-1 1];
% R5: RL -> R
N([2 1],5)=[-1 1];
% R6: G -> Ga + Gbg 
N([3 4 5],6)= [-1 1 1];
% R7: Ga -> Gd
N([4 6],7)=[-1 1];
% R8: Gd + Gbg -> G
N([6 5 3],8)=[-1 -1 1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
umod.N = sparse(N);

Ncells = umod.Ncells; 
Ndofs=Ncells*Mspecies;

% Use the subdomain vector sd to find those dofs (degress of freedom)
% that are on the membrane (pm) and in the cytosol (cyt)  
pm  = find(umod.sd == 2);
cyt = find(umod.sd == 1);

% Setup Initial Conditions
u0 = zeros(Mspecies,Ncells);
u0(1,pm)=round((10000/length(pm)));
u0(3,pm)=round((10000/length(pm)));

umod.u0 = u0;

% Setup  sample times (every 10 seconds)
umod.tspan = 0:10:100;



%evaluate the gradient expression at the meshpoints
if ~isjava(umod.comsol)
    data = postinterp(umod.comsol,'Gradient',umod.comsol.mesh.p);
    dofs = xmeshinfo(umod.comsol,'out','dofs');
    data = data(dofs.nodes(1:Mspecies:end));
else
    xmi = mphxmeshinfo(umod.comsol);
    coords = xmi.nodes.coords;
    data = mphinterp(umod.comsol,{'Gradient'}, 'coord', coords);
    dofs = xmi.dofs;
    data = data(dofs.nodes(1:Mspecies:end)+1);
    data = data';
end
% permutation to correct order for urdme

umod.data = data;

% Verbose mode on solver. 
%umod.options = [2 NaN 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We need to modify the diffusion matrix such that the membrane
% bound species only diffuse on the membrane. We achieve this by
% zeroing out the elements corresponding to diffusion in the 
% cytosol. 

% For all species, flag all dofs in the cytosol for removal. 
ixremove = [];
for s = 1:Mspecies
  ixremove = [ixremove; Mspecies*(cyt-1)+s];
end

% Decompose the sparse matrix. 
[i,j,s] = find(umod.D');

% And set all elements in the diffusion matrix corresponding to 
% to the cytosol to zero.

ixremove = [find(ismember(i,ixremove)); find(ismember(j,ixremove))];
i(ixremove) = [];
j(ixremove) = [];
s(ixremove) = [];

% Reassemble the sparse matrix and adjust the diagonal entries. 
ixkeep = find(s > 0);
D = sparse(i(ixkeep),j(ixkeep),s(ixkeep),Ndofs,Ndofs);
d = full(sum(D,2));
D = D+sparse(1:Ndofs,1:Ndofs,-d);

umod.D = D';
