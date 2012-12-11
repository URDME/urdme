function umod = urdme2comsol(umod,U,tspan,verbose)
%URDME2COMSOL Update comsol model in umod.comsol with result matrix U
%   UMOD = URDME2COMSOL(UMOD,U,VERBOSE) Inserts solution matrix U into
%   Umod.comsol
%
%   --INPUT--
%   UMOD - Structure containing URDME simulation parameters and the
%       comsol model in umod.comsol.
%   U - Solution matrix produced by URDME core. U are in absolute integers
%       (number of species), not concentration. URDME2COMSOL turns values to 
%       consentration before saving to the model object. 
%   tspan - Times span vector that matches the size of U. Defaults to umod.tspan.  
%   VERBOSE - Silent execution
%
%   --OUTPUT--
%   UMOD - Same as input but with U added to Solution 'sol1' of umod.cosmol
%
%   See also COMSOL2URDME, URDME.

% P.Bauer       2012-08-30
% V.Gerdin      2012-02-27 (rdme2mod)
% A. Hellander  2010-05-04 (rdme2fem)
% J. Cullhed    2008-08-06 (rdme2fem)

if(isfield(umod,'comsol'))
    is4x = isjava(umod.comsol);
else
    is4x = 0;
end

if nargin < 3
    tspan = umod.tspan;
end

if nargin < 4
    verbose = 0;
end
  

if is4x
  %extend mesh from model
  xmi = mphxmeshinfo(umod.comsol);

  % RESHAPE U TO FIT MODEL

  % fundamentals
  nodes = xmi.nodes;
  [Mspecies,Ncells] = size(nodes.dofs);
  Ndofs = Mspecies*Ncells;

  % permutation of species
  if isfield(umod,'species')
      fieldnames = umod.species;
  else
      fieldnames = xmi.fieldnames;
  end

  % Fixed for Comsol 4.3
  modName = char(umod.comsol.modelNode.tags());
  if isempty(strfind(char(fieldnames(1)),modName))
    fieldnames=strcat(modName,'.',fieldnames);
  end
  [foo,sp] = ismember(fieldnames,nodes.dofnames);
  isp(sp) = 1:numel(sp); %invert order

  % Scale U
  U = reshape(U,Mspecies,Ncells,[]);

  % Change from absolute figures to concentration
  % suboptimal, but should work ok
  vol = umod.vol;

  V = repmat(vol(:)'*6.022e23,[Mspecies 1]);
  if ~exist(tspan)
    tspan = umod.tspan;
  end
  for i = 1:numel(tspan)
      U(:,:,i) = U(:,:,i)./V;
  end

  % Finally reshape U in order of Dofs
  U = reshape(U(isp,:,:),Ndofs,[]);

  % Store U in model
  if verbose>0
    disp('Saving U to Model')
  end

  for i=1:1:length(tspan)
      umod.comsol.sol('sol1').setU(i,U(:,i));
      if verbose>1
        disp(sprintf('Saved %2.0f%%',i*100/length(tspan)))
      end
  end

  umod.comsol.sol('sol1').setPVals(tspan);
  umod.comsol.sol('sol1').createSolution;

else

  % extend mesh
  xmesh = meshextend(umod.comsol);

  % fundamentals
  nodes = xmeshinfo(xmesh,'out','nodes');
  [Mspecies,Ncells] = size(nodes.dofs);
  Ndofs = Mspecies*Ncells;

  % permutation of species
  [foo,sp] = ismember(umod.comsol.equ.dim,nodes.names);
  isp(sp) = 1:numel(sp);

  % scale and permute back
  U = reshape(U,Mspecies,Ncells,[]);

  % suboptimal, but should work ok:
  vol = umod.vol;
  V = repmat(vol(:)'*6.022e23,[Mspecies 1]);
  for i = 1:numel(tspan)
    U(:,:,i) = U(:,:,i)./V;
  end

  % final shape
  U = reshape(U(isp,:,:),Ndofs,[]);
  
  % Store U in model
  if verbose>0
    disp('Saving U to Model')
  end

  % final call
  umod.comsol.sol = femsol(U, 'tlist', tspan);  
  
end

if verbose>0
  disp('Saving: Done')
end