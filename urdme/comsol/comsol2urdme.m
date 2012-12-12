function umod = comsol2urdme(fem,varargin)
%COMSOL2URDME Create an URDME-structure from a COMSOL multiphysics model
%   UMOD = COMSOL2URDME(FEM,VARARGIN)
%
%   --INPUT--
%   FEM - Comsol Multiphysics 3.5a or 4.x model file
%   VERBOSE - Silent execution
%
%   --OUTPUT--
%   UMOD - Urdme 1.2 structure with fem-model added to umod.comsol
%
%   Generalized subdomains (non-empty argument SDUSED). 
%   The generalized subdomain vector is given a default value of one. Other 
%   values can be defined by a special expression (rdme_sd) which should be 
%   defined on every subdomain, boundary and point. The global expression 
%   rdme_sdlevel should be defined as the lowest dimension that rdme_sd is 
%   defined for.
%

% P. Bauer     2012-08-30
% V. Gerdin    2012-02-08 (mod2rdme)
% S. Engblom   2008-08-01 (fem2rdme)
% A. Hellander 2009-11-18 (fem2rdme)
is4x = isjava(fem);
if nargin > 2
    verbose = varargin{1};
else
    verbose = 0;
end

if is4x  
  if verbose>0    
    fprintf('Identified fem-model as Comsol 4.x.\n');
  end
      
  % Get extended mesh data
  if verbose>0
    disp('Start mphxmeshinfo')
  end

  xmi = mphxmeshinfo(fem); 
  if verbose>0
    disp('Done')
  end

  % Get number of species and number of cells.
  [Mspecies,Ncells] = size(xmi.nodes.dofs);
  Ndofs = Ncells*Mspecies;
  if Ndofs ~= xmi.ndofs
      warning(['Number of degrees of freedoms mismatch. ' ...
               'Check the boundary conditions.']);
  end
  % permutation of species
  [foo,sp] = ismember(xmi.fieldnames,xmi.nodes.dofnames);

  % permutation fem -> rdme
  f2r = reshape(1:Ndofs,Mspecies,Ncells);
  f2r = reshape(f2r(sp,:),[ ],1);

  % Create tspan & u0 placeholder, must be redefined in model.m
  tspan = zeros(1,0);
  u0 = zeros(1,0);

  % Assemble matrices
  if verbose>0
    disp('Start mphmatrix')
  end
  DM = mphmatrix(fem,'sol1','out',{'K','D'});
  if verbose>0
    disp('mphmatrix: Done')
  end
  D = DM.K; % URDME: D "Diffusion matrix"  Comsol: K "Stiffness matrix"
  M = DM.D; % URDME: M "Mass matrix"  Comsol: D "Damping matrix"

  % The (lumped) mass matrix gives the element volume
  vol = full(sum(M,2));

  % explicitly invert the lumped mass matrix and filter the diffusion matrix
  [i,j,s] = find(D);
  diagsum1 = sum(diag(D)./vol);
  s = -s./vol(i);
  ixkeep = find(s > 0);
  i = i(ixkeep); j = j(ixkeep); s = s(ixkeep);

  % permute and rebuild
  if2r(f2r) = 1:Ndofs;
  [i,j,p] = find(sparse(if2r(i),if2r(j),1:numel(s)));
  D = sparse(i,j,s(p),Ndofs,Ndofs);
  d = full(sum(D,2));
  D = D+sparse(1:Ndofs,1:Ndofs,-d);
  diagsum2 = sum(d);

  % check if the difference is too big
  if abs(diagsum2-diagsum1) > abs(diagsum1)*0.10
    warning('Many off diagonal negative elements in D.');
  end
  % need only volume of unique cells
  vol = vol(f2r);
  vol = vol(1:Mspecies:end);

  % Generalized subdomains. To support more than one subdomain in 
  % the models, the global expression "rdme_sdlevel" has to 
  % be defined.

  sd = ones(1,Ncells);

  % Find rdme_sdlevel value if present
  rdme_sdlvl = -1;
  if verbose>0
    disp('Start mphglobal')
  end
  try
      rdme_sdlvl = mphglobal(fem,'rdme_sdlevel','solnum','end');
  catch
      warning('Global variable "rdme_sdlevel" not set. No subdomains will be used.');
  end
  if verbose>0
    disp('Done')
  end

  if rdme_sdlvl > -1 % i.e. if we found a valid rdme_sdlvl
    % initial sd values, set all to 0
    sd = zeros(1,Ncells);

    % sdim - the maximum number of dimensions to run simulation
    sdim = length(xmi.elements.meshtypes)-1;
    % edim2str - translation from numeric edim to string
    edim2str = {'vtx' 'edg' 'tri' 'tet'};

    for edim = sdim:-1:rdme_sdlvl

      % enode_idx - the indices of the nodes involved for a given edim
      enode_idx = unique(xmi.elements.(edim2str{edim+1}).nodes+1);
      % add 1 to every index as 4.2 indices start from 0

      % coords - x,y,z coordinates of each involved node 
      coords = xmi.nodes.coords(:,enode_idx);
      if verbose>0
        disp('Start mphinterp')
      end
      try        
        rdme_sd = mphinterp(fem, 'rdme_sd', 'coord', coords, 'edim', edim, 'solnum', 1);
      catch
          error(['To use more than une subdomain, you need to specify' ...
                 'the expression rdme_sd, see manual.'])
      end
      if verbose>0
        disp('Done')
      end
      % Set sd for involved nodes
      sd(enode_idx) = rdme_sd;

    end
    % sd is currently ordered according to comsol nodes. We need to permute
    %     it to be ordered according to comsol dofs (like everything else).
    sd = sd(xmi.dofs.nodes(1:Mspecies:end)+1);
    % add 1 to every index as 4.2 indices start from 0
  end

  % as of the new Comsol API the species' names are refered to as
  % "modelID.speciesName", ex "mod1.MinD_c_atp". The following code
  % removes the new prefix
  tags = fem.modelNode.tags;
  species = xmi.fieldnames;

  for i = [1:size(xmi.fieldnames)]
      if strfind(xmi.fieldnames{i},[char(tags(1)) '.']) == 1
          fullName = char(xmi.fieldnames{i});
          species{i} = fullName(size(char(tags(1)),2)+2:end);
      end
  end
      
      
else
   if verbose>0    
    fprintf('Identified fem-model as Comsol 3.5a.\n');
  end
      
  % extend mesh
  if verbose>0
    disp('Start meshextend')
  end
  try
  fem.xmesh = meshextend(fem);
  catch
      umod=fem;
      return;
  end
  if verbose>0
    disp('Done')
  end

  % fundamentals
  if verbose>0
    disp('Start xmeshinfo')
  end
  [nodes,dofs] = xmeshinfo(fem.xmesh,'out',{'nodes' 'dofs'});
  if verbose>0
    disp('Done')
  end
  [Mspecies,Ncells] = size(nodes.dofs);
  Ndofs = Ncells*Mspecies;
  if Ndofs ~= flngdof(fem)
    warning(['Number of degrees of freedoms mismatch. ' ...
            'Check the boundary conditions.']);
  end
  % permutation of species
  [foo,sp] = ismember(fem.equ.dim,nodes.names);

  % permutation fem -> rdme
  f2r = reshape(1:Ndofs,Mspecies,Ncells);
  f2r = reshape(f2r(sp,:),[ ],1);

  % create tspan (if any)
  tspan = zeros(1,0);
  if isfield(fem,'sol') && isfield(fem.sol,'tlist')
    tspan = fem.sol.tlist;
  end

  % u0
  if verbose>0
    disp('Start asseminit')
  end
  u0 = asseminit(fem);
  if verbose>0
    disp('Done')
  end
  u0 = reshape(u0.u(f2r),Mspecies,Ncells);

  % assemble matrices
  if verbose>0
    disp('Start assemble')
  end
  [D,M] = assemble(fem,'out',{'K','D'});
  if verbose>0
    disp('Done')
  end

  % the (lumped) mass matrix gives the element volume
  vol = full(sum(M,2));

  % explicitly invert the lumped mass matrix and filter the diffusion matrix
  [i,j,s] = find(D);
  diagsum1 = sum(diag(D)./vol);
  s = -s./vol(i);
  ixkeep = find(s > 0);
  i = i(ixkeep); j = j(ixkeep); s = s(ixkeep);

  % permute and rebuild
  if2r(f2r) = 1:Ndofs;
  [i,j,p] = find(sparse(if2r(i),if2r(j),1:numel(s)));
  D = sparse(i,j,s(p),Ndofs,Ndofs);
  d = full(sum(D,2));
  D = D+sparse(1:Ndofs,1:Ndofs,-d);
  diagsum2 = sum(d);

  % check if the difference is too big
  if abs(diagsum2-diagsum1) > abs(diagsum1)*0.10
    warning('Many off diagonal negative elements in D.');
  end

  % switch to absolute numbers
  u0real = 6.022e23*vol'*u0(:);
  u0 = round(6.022e23*reshape(vol,Mspecies,Ncells).*u0);

  % check if the difference is too big
  if abs(u0real-sum(u0(:))) > u0real*1e-2
    warning('The initial concentration is not accurately represented.');
  end

  % need only volume of unique cells
  vol = vol(f2r);
  vol = vol(1:Mspecies:end);

  % Generalized subdomains. To support more than one subdomain in 
  % the models, the global expression "rdme_sdlevel" has to 
  % be defined.

  sd = ones(1,Ncells);
  if isfield(fem,'globalexpr')
      i = ismember(fem.globalexpr(1:2:end),'rdme_sdlevel');
      [sdused,j]=max(i);
      if sdused
      temp=fem.globalexpr(2*j);
      rdme_sdlevel = str2num(temp{1});
      sd = zeros(1,Ncells);
      if verbose>0
        disp('Start xmeshinfo')
      end
      els = xmeshinfo(fem,'out','elements');
      if verbose>0
        disp('Done')
        disp('Start xmeshinfo')
      end
      sdim = max(xmeshinfo(fem,'out','edim'));
      if verbose>0
        disp('Done')
      end
      for edim = sdim:-1:rdme_sdlevel
          enodes = els{edim+1}.nodes;
          if edim == 0
              spoint = [0];
          else
              spoint = [zeros(edim,1) eye(edim)];
          end
          try
              if verbose>0
                disp('Start posteval')
              end
              rdme_sd = posteval(fem,'rdme_sd','spoint',spoint,'edim',edim);
              if verbose>0
                disp('Done')
              end
          catch
              error(['To use more than une subdomain, you need to specify' ...
              'the expression rdme_sd, see manual.'])
          end
          sd(enodes(:,rdme_sd.elind(1:edim+1:end))) = rdme_sd.d;
      end
      sd = sd(dofs.nodes(1:Mspecies:end));
      end
  end
end  

% Create URDME 1.2 structure
umod = struct();
umod.version = '1.2';
if(isjava(fem))
  umod.mph     = '4.x';
else
  umod.mph     = '3.5a';	
end

umod.tspan   = tspan;
umod.u0      = u0;
umod.D       = D';
umod.vol     = vol;
umod.sd      = sd;
umod.Ncells  = Ncells;
umod.sopts   = [];
umod.comsol  = fem;

	