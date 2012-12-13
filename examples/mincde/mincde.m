% Model file for the MINCDE example.
%
% P. Bauer 2012-11-22
% A. Hellander 2010-06-07
%
%

function umod = mincde(umod)

Mspecies   = 5;
Mreactions = 5;
Ncells = umod.Ncells; 
Ndofs = Mspecies*Ncells;

% Stoichiometry matrix. Every column corresponds to a reaction.
umod.N=sparse([-1  -1  0  0  1 ;...
                1   1 -1  0  0 ;...
                0   0 -1  1  0 ;...
                0   0  1 -1  0 ;...
                0   0  0  1 -1]);
                
% Dependency graph. The first Mspecies columns tells which propensities
% needs to be updated that species diffuses. The following Mreactions
% columns does the same thing but for reaction events.
umod.G=sparse([1 0 0 0 0 1 1 0 0 1;...
               1 1 0 0 0 1 1 1 0 1;...
               0 1 1 0 0 1 1 1 1 0;...
               0 0 0 1 0 0 0 1 1 0;...
               0 0 0 0 1 0 0 0 1 1;]);
                

u0 = zeros(Mspecies,Ncells);

% Set "random" initial distribution. This can be achieved in several
% ways, below is an example of one way of doing it. 
nMinD = 4002;
nMinE = 1040;

ind = floor(Ncells*rand(1,nMinD))+1;
u0(5,:) = full(sparse(1,ind,1,1,Ncells));

ind = floor(Ncells*rand(1,nMinE))+1;
u0(3,:) = full(sparse(1,ind,1,1,Ncells));

umod.u0 = u0;

% The time span. The state of the trajectory will be stored at all
% time points in tspan. 
umod.tspan = 0:400;

% We need to modify the diffusion matrix such that the membrane
% bound species only diffuse on the membrane. We achieve this by
% zeroing out the elements corresponding to diffusion in the 
% cytosol. 

% Use the subdomain vector sd to find the dofs (degress of freedom)
% that are on the membrane (pm) and in the cytosol (cyt)  
pm  = find(umod.sd == 2);
cyt = find(umod.sd == 1);

% For MinD_m (2) and MinDE (4), flag all dofs in the cytosol for 
% removal. 
ixremove = [];
for s = [2 4]
  ixremove = [ixremove; Mspecies*(cyt-1)+s];
end


D = umod.D';

% Decompose the sparse matrix. 
[i,j,s] = find(D);

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

% One of the propensity functions needs to be scaled by 1/h, where
% h is the local lengthscale of the mesh (for a Cartesian mesh, 
% this corresponds to the side lenghts of the cubes). We can obtain that
% information by prompting Comsol for that data. Using the data vector 
% (input to the propensity functions) we pass h to each propesity function.

% "postinterp" is a built-in Comsol function that evaluates any (valid)
% expression at a set of specified points. The predefined 
% expression 'h' can be used to obtain the local length of each subvolume. 

if iscmp4x(umod.comsol) %Comsol 4.x
  xmi = mphxmeshinfo(umod.comsol);
  umod.data = mphinterp(umod.comsol,'h','coord', xmi.dofs.coords(:,1:Mspecies:end), 'solnum', 1);
  umod.data = umod.data(xmi.dofs.nodes(1:Mspecies:end)+1);
else
  dofs = xmeshinfo(umod.comsol,'Out','dofs');
  umod.data = postinterp(umod.comsol,'h',dofs.coords(:,1:Mspecies:end));
  umod.data = umod.data(dofs.nodes(1:Mspecies:end));
end



