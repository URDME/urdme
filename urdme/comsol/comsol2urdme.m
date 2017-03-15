function umod = comsol2urdme(fem,verbose)
%COMSOL2URDME Create URDME-structure from a COMSOL multiphysics model.
%   UMOD = COMSOL2URDME(COMSOL,VERBOSE) creates an URDME-structure
%   UMODE from the Comsol Java object COMSOL under the verbosity
%   VERBOSE.
%
%   See also URDME, URDME2COMSOL.

% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% P. Bauer 2012-08-30 (Revision)
% V. Gerdin 2012-02-08 (mod2rdme)
% S. Engblom 2008-08-01 (fem2rdme)

% Generalized subdomains (non-empty argument SDUSED). The generalized
% subdomain vector is given a default value of one. Other values can
% be defined by a special expression (rdme_sd) which should be defined
% on every subdomain, boundary and point. The global expression
% rdme_sdlevel should be defined as the lowest dimension that rdme_sd
% is defined for.

if nargin < 2, verbose = 0; end

% Get extended mesh data
l_info(verbose,1,'Start mphxmeshinfo...');
xmi = mphxmeshinfo(fem); 
l_info(verbose,1,' ...done.\n');

% Get number of species and number of cells.
[Mspecies,Ncells] = size(xmi.nodes.dofs);
Ndofs = Ncells*Mspecies;
if Ndofs ~= xmi.ndofs
  warning(['Number of degrees of freedoms mismatch. ' ...
           'Check the boundary conditions.']);
end
% permutation of species
[foo,sp] = ismember(xmi.fieldnames,xmi.nodes.dofnames);

% permutation comsol -> urdme
f2r = reshape(1:Ndofs,Mspecies,Ncells);
f2r = reshape(f2r(sp,:),[ ],1);

% Assemble matrices
l_info(verbose,1,'Start mphmatrix...');
DM = mphmatrix(fem,'sol1','out',{'K','D'});
l_info(verbose,1,' ...done.\n');
D = DM.K; % URDME: D "Diffusion matrix" Comsol: K "Stiffness matrix"
M = DM.D; % URDME: M "Mass matrix"      Comsol: D "Damping matrix"

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

% Generalized subdomains. To support more than one subdomain in the
% models, the global expression "urdme_sdlevel" has to be defined.
sd = ones(1,Ncells);

% Find rdme_sdlevel value if present
urdme_sdlvl = -1;
l_info(verbose,1,'Start mphglobal...')
try
  urdme_sdlvl = mphglobal(fem,'urdme_sdlevel','solnum','end');
catch
  l_info(verbose,1,['Global variable "urdme_sdlevel" not set. ' ...
                    'No subdomains will be used.']);
end
l_info(verbose,1,' ...done.\n')

% if we found a valid urdme_sdlvl initial sd values, set all to 0
if urdme_sdlvl > -1
  sd = zeros(1,Ncells);

  % sdim - the maximum number of dimensions to run simulation
  sdim = length(xmi.elements.meshtypes)-1;
  % edim2str - translation from numeric edim to string
  edim2str = {'vtx' 'edg' 'tri' 'tet'};
  for edim = sdim:-1:urdme_sdlvl
    % enode_idx - the indices of the nodes involved for a given edim
    enode_idx = unique(xmi.elements.(edim2str{edim+1}).nodes+1);
    % add 1 to every index as 4.2 indices start from 0

    % coords - x,y,z coordinates of each involved node 
    coords = xmi.nodes.coords(:,enode_idx);
    l_info(verbose,1,'Start mphinterp...');
    try        
      urdme_sd = mphinterp(fem,'urdme_sd','coord',coords,'edim',edim,'solnum',1);
    catch
      error(['To use more than une subdomain, you need to specify ' ...
             'the expression urdme_sd, see manual.'])
    end
    l_info(verbose,1,' ...done.\n');
    % Set sd for involved nodes
    sd(enode_idx) = urdme_sd;
  end
  % sd is currently ordered according to comsol nodes. We need to
  % permute it to be ordered according to comsol dofs (like everything
  % else).
  sd = sd(xmi.dofs.nodes(1:Mspecies:end)+1);
  % add 1 to every index as 4.2 indices start from 0
end

% In the Comsol API the species' names are refered to as
% "modelID.speciesName", e.g., "mod1.MinDcytATP". The following code
% removes the prefix.
tags = fem.modelNode.tags;
species = xmi.fieldnames;

for i = 1:size(xmi.fieldnames)
  if strfind(xmi.fieldnames{i},[char(tags(1)) '.']) == 1
    fullName = char(xmi.fieldnames{i});
    species{i} = fullName(size(char(tags(1)),2)+2:end);
  end
end

% finalize the Comsol-part of the URDME structure
umod = struct('D',D', ...
              'vol',vol, ...
              'sd',sd, ...
              'comsol',fem);

%-------------------------------------------------------------------------
function l_info(level,lim,msg)
%L_INFO Display information.
%   L_INFO(LEVEL,LIM,MSG) Displays the message MSG whenever LEVEL >=
%   LIM.

if level >= lim, fprintf(msg); end

%-------------------------------------------------------------------------
