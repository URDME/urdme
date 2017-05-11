function umod = comsol2urdme(fem,verbose)
%COMSOL2URDME Create URDME-structure from a COMSOL multiphysics model.
%   UMOD = COMSOL2URDME(COMSOL,VERBOSE) creates an URDME-structure
%   UMODE from the Comsol Java object COMSOL under the verbosity
%   VERBOSE.
%
%   See also URDME, URDME2COMSOL.

% S. Engblom 2017-05-08 (Revision, conceptual change in assembly)
% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% P. Bauer 2012-08-30 (Revision)
% V. Gerdin 2012-02-08 (mod2rdme)
% S. Engblom 2008-08-01 (fem2rdme)

% Generalized subdomains (non-empty argument SDUSED). The generalized
% subdomain vector is given a default value of one. Other values can
% be defined by a special expression (urdme_sd) which should be
% defined on every subdomain, boundary, and point of the geometry. The
% global expression urdme_sdlevel should be defined as the lowest
% dimension that urdme_sd is defined for.

if nargin < 2, verbose = 0; end

% get extended mesh data
l_info(verbose,1,'Start mphxmeshinfo...');
try
  xmi = mphxmeshinfo(fem);
catch ME
  % Comsol-problems: sometimes the solution is lost and a model re-run
  % is the most immediate solution
  if strcmp(ME.message,'The model does not contain any solutions')
    warning(['No solution found in the Comsol object. ' ...
             'Will now try to execute model.sol.run().']);
    fem.sol.run()
    xmi = mphxmeshinfo(fem);
  end
end
l_info(verbose,1,' ...done.\n');

% number of species and number of cells
[Mspecies,Ncells] = size(xmi.nodes.dofs);
Ndofs = Ncells*Mspecies;
if Ndofs ~= xmi.ndofs
  warning(['Number of degrees of freedoms mismatch. ' ...
           'Check the boundary conditions.']);
end

% permutation of species, comsol --> urdme
[foo,sp] = ismember(xmi.fieldnames,xmi.nodes.dofnames);
f2r = reshape(1:Ndofs,Mspecies,Ncells);
f2r = reshape(f2r(sp,:),[ ],1);

% assemble matrices
l_info(verbose,1,'Start mphmatrix...');
DM = mphmatrix(fem,'sol1','out',{'K','D'});
l_info(verbose,1,' ...done.\n');
D = DM.K; % URDME: D "Diffusion matrix" Comsol: K "Stiffness matrix"
M = DM.D; % URDME: M "Mass matrix"      Comsol: D "Damping matrix"

% the (lumped) mass matrix gives the element volume
vol = full(sum(M,1))';
vol = vol(f2r);

% multiply from the right with the inverse lumped mass matrix and
% filter the diffusion matrix
[i,j,s] = find(D(f2r,f2r));
diagsum1 = sum(s(i == j)./vol);
s = -s./vol(j); % from the right (i.e. columns) + sign convention
ixkeep = find(s > 0);
i = i(ixkeep); j = j(ixkeep); s = s(ixkeep);

% rebuild
D = sparse(i,j,s,Ndofs,Ndofs);
d = full(sum(D,1));
D = D+sparse(1:Ndofs,1:Ndofs,-d);
diagsum2 = sum(d);

% check if the sum difference is too big
if abs(diagsum2-diagsum1) > abs(diagsum1)*0.10
  warning('Many off diagonal negative elements in D.');
end

% save only the volume of unique cells
vol = vol(1:Mspecies:end);

% generalized subdomains: to support more than one subdomain in the
% models, the global expression "urdme_sdlevel" has to be defined
sd = ones(1,Ncells);

% find urdme_sdlevel value if present
urdme_sdlvl = -1;
l_info(verbose,1,'Start mphglobal...')
try
  urdme_sdlvl = mphglobal(fem,'urdme_sdlevel','solnum','end');
catch
  l_info(verbose,1,['Global variable "urdme_sdlevel" not set. ' ...
                    'No subdomains will be used.']);
end
l_info(verbose,1,' ...done.\n')

% if we found a valid urdme_sdlvl set initial sd values to 0
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
      urdme_sd = mphinterp(fem,'urdme_sd', ...
                           'coord',coords,'edim',edim,'solnum',1);
    catch
      error(['To use more than one subdomain, you need to specify ' ...
             'the expression urdme_sd, see manual.'])
    end
    l_info(verbose,1,' ...done.\n');
    % Set sd for involved nodes
    sd(enode_idx) = urdme_sd;
  end
  % sd is currently sized and ordered according to Comsol dofs and so we
  % need to select only unique voxels
  sd = sd(xmi.dofs.nodes(1:Mspecies:end)+1);
  % add 1 to every index as 4.2 indices start from 0
end

% not currently used
% $$$ % in the Comsol API the species' names are refered to as
% $$$ % "modelID.speciesName", e.g., "mod1.MinDcytATP"
% $$$ tags = fem.modelNode.tags;
% $$$ species = xmi.fieldnames;
% $$$ for i = 1:size(xmi.fieldnames)
% $$$   if strfind(xmi.fieldnames{i},[char(tags(1)) '.']) == 1
% $$$     fullName = char(xmi.fieldnames{i});
% $$$     species{i} = fullName(size(char(tags(1)),2)+2:end);
% $$$   end
% $$$ end

% finalize the Comsol-part of the URDME structure
umod = struct('D',D, ...
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
