function umod = urdme2pde(umod)
%URDME2PDE Update .pde-field with URDME result.
%   UMOD = URDME2PDE(UMOD) inserts the trajectory (UMOD.U,UMOD.tspan)
%   into UMOD.pde under the verbosity VERBOSE.
%
%   See also URDME, PDE2URDME, COMSOL2URDME.

% S. Engblom 2017-02-21

% fundamentals
Mspecies = size(umod.N,1);
Ncells = numel(umod.vol);

% change from absolute figures to concentration
U = reshape(umod.U,Mspecies,Ncells,[]);
V = repmat(umod.vol(:)'*6.022e23,[Mspecies 1]);
for i = 1:numel(umod.tspan)
  U(:,:,i) = U(:,:,i)./V;
end

% store U in model.pde.U
umod.pde.U = U;
