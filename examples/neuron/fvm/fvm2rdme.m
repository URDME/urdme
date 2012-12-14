% Assemble first order finite volume discretization for a convection term. 

A = function fvm2rdme(fem,file)

% extend mesh
if ~isfield(fem,'xmesh')
  fem.xmesh = meshextend(fem);
end

% fundamentals
[nodes,dofs] = xmeshinfo(fem.xmesh,'out',{'nodes' 'dofs'});
[Mspecies,Ncells] = size(nodes.dofs);
Ndofs = Ncells*Mspecies;

if Ndofs ~= flngdof(fem)
  warning(['Number of degrees of freedoms mismatch. ' ...
           'Check the boundary conditions.']);
end



