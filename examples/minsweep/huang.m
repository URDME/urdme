function umod = huang(umod)
%HUANG Matlab model file for the MinD/MinE-sweep example.
%   The model is taken from Huang et al. (2003).
%
%   MinD_c_atp -> MinD_m
%   MinD_c_atp + MinD_m -> 2MinD_m (cooperative binding)
%   MinD_m + MinD_e -> MinDE
%   MinDE -> MinD_c_adp + MinD_e
%   MinD_c_adp -> MinD_c_atp
%   MinDE + MinD_c_atp -> MinD_m + MinDE
%
%   Ordering of species:
%     [MinD_c_atp MinD_m MinD_e MinDE MinD_c_adp].

% V. Gerdin 2012-03-09 (Urdme 4.2 support)
% S. Engblom 2011-06-01 (Revision)
% A. Hellander 2010-06-07

% dimensions
Mspecies = 5;
Mreactions = 6;
Ncells = umod.Ncells; 
Ndofs = Mspecies*Ncells;

% stoichiometric matrix
umod.N = sparse([-1 -1  0  0  1 -1; ...
                       1  1 -1  0  0  1; ...
                       0  0 -1  1  0  0; ...
                       0  0  1 -1  0  0; ...
                       0  0  0  1 -1  0]);

% dependency graph Mreactions-by-(Mreactions+Mspecies)
umod.G = sparse([1 0 0 0 0 1 1 0 0 1 1; ...
                      1 1 0 0 0 1 1 1 0 1 1; ...
                      0 1 1 0 0 1 1 1 1 0 1; ...
                      0 0 0 1 0 0 0 1 1 0 0; ...
                      0 0 0 0 1 0 0 0 1 1 0; ...
                      1 1 1 1 1 1 1 1 1 1 1]);

% build a random initial distribution
nMinD_e = 1575;
nMinD_c_adp = 4500;
u0 = zeros(Mspecies,Ncells);
u0(3,:) = full(sparse(1,ceil(Ncells*rand(1,nMinD_e)),3,1,Ncells));
u0(5,:) = full(sparse(1,ceil(Ncells*rand(1,nMinD_c_adp)),1,1,Ncells));
umod.u0 = u0;

% default sampling interval
umod.tspan = 0:100:2000;

% assigne the value of 'h' (local element size) to the data vector
if iscmp4x(umod.comsol) % Comsol 4.x
    xmi = mphxmeshinfo(umod.comsol);
    umod.data = mphinterp(umod.comsol,'h','coord',xmi.dofs.coords(:,1:Mspecies:end));
    umod.data = umod.data(xmi.dofs.nodes(1:Mspecies:end)+1);
    umod.data = umod.data';  
else
    dofs = xmeshinfo(umod.comsol,'Out','dofs');
    umod.data = postinterp(umod.comsol,'h',umod.comsol.mesh.p); 
    umod.data = umod.data(dofs.nodes(1:Mspecies:end));
end

