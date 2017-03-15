function umod = mincde(umod)
% Model file for the MINCDE example.

% S. Engblom 2017-02-20 (Revision, URDME 1.3, Comsol 5.2)
% P. Bauer 2012-11-22
% A. Hellander 2010-06-07
% J. Cullhed 2008-08-14

% create the reaction network
r1 = 'MinDcytATP > sd == 2 ? kd*MinDcytATP/ldata[0] : 0.0 > MinDmem';
r2 = ['MinDcytATP + MinDmem > kdD*MinDcytATP*MinDmem/(1000.0*NA*vol) > ' ...
      'MinDmem+MinDmem'];
r3 = 'MinE + MinDmem > kde*MinE*MinDmem/(1000.0*NA*vol) > MinDE';
r4 = 'MinDE > ke*MinDE > MinDcytADP + MinE';
r5 = 'MinDcytADP > kp*MinDcytADP > MinDcytATP';

species = {'MinDcytATP' 'MinDmem' 'MinE' 'MinDE' 'MinDcytADP'};
rates = {'NA' 6.022e23 'kd' 1.25e-8 'kdD' 9.0e6 'kde' 5.56e7 ...
         'ke' 0.7 'kp' 0.5};
[~,umod.N,umod.G] = rparse({r1 r2 r3 r4 r5},species,rates,'fange.c');

Mspecies = size(umod.N,1);
Ncells = numel(umod.sd);
Ndofs = Mspecies*Ncells;

% Set "random" initial distribution. This can be achieved in several
% ways, below is one example of doing it.
u0 = zeros(Mspecies,Ncells);
nMinE = 1040;
nMinD = 4002;

ind = floor(Ncells*rand(1,nMinE))+1;
u0(3,:) = full(sparse(1,ind,1,1,Ncells));

ind = floor(Ncells*rand(1,nMinD))+1;
u0(5,:) = full(sparse(1,ind,1,1,Ncells));

umod.u0 = u0;

% The state of the trajectory will be stored at all time points in
% tspan.
umod.tspan = 0:200;

% We need to modify the diffusion matrix such that the membrane bound
% species only diffuse on the membrane. We achieve this by zeroing out
% the elements corresponding to diffusion in the cytosol.

% Use the subdomain vector sd to find the dofs (degress of freedom)
% that are on the membrane (pm) and in the cytosol (cyt)
cyt = find(umod.sd == 1);
pm  = find(umod.sd == 2);

% For MinDmem (#2) and MinDE (#4), flag all dofs in the cytosol for
% removal.
ixremove = [];
for s = [2 4]
  ixremove = [ixremove; Mspecies*(cyt-1)+s];
end

% Decompose the sparse matrix. 
D = umod.D';
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

% One of the propensity functions needs to be scaled by 1/h, where h
% is the local lengthscale of the mesh (for a Cartesian mesh, this
% corresponds to the side lenghts of the cubes). We can obtain that
% information by prompting Comsol for that data. Using the ldata
% vector (input to the propensity functions) we pass h to each
% propensity function.
%
% "mphinterp" is a built-in Comsol function that evaluates an
% expression at a set of specified points. The predefined expression
% 'h' can be used to obtain the local length of each subvolume.
xmi = mphxmeshinfo(umod.comsol);
umod.ldata = mphinterp(umod.comsol,'h','coord', ...
                       xmi.dofs.coords(:,1:Mspecies:end),'solnum',1);
umod.ldata = umod.ldata(xmi.dofs.nodes(1:Mspecies:end)+1);
