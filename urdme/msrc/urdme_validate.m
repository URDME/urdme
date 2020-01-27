function urdme_validate(umod)
% URDME_VALIDATE Validate URDME struct.
%    URDME_VALIDATE(UMOD) attempts to validate the fields in UMOD and
%    throws an error if a field is missing or has an incorrect format.
%
%    See also URDME.

% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2019-11-12 (Revision, multiple seeds)
% S. Engblom 2018-02-10 (Revision, Nreplicas syntax)
% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% J. Cullhed 2008-06-18 (rdme.m)

% fields checked
req = {'tspan' 'u0' 'D' 'N' 'G' 'vol' 'sd' ...
      'report' 'seed' 'inline_propensities' ...
      'ldata' 'gdata' 'solverargs' 'makeargs'};
% check that fields exist first
for i = 1:numel(req)
  if ~isfield(umod,req{i})
    error(['No field .' req{i} '.']);
  end
end

% Check tspan.
tspan = reshape(umod.tspan,1,[]);
if issparse(tspan) || ~isa(tspan,'double')
  error('Field .tspan must be a double vector.');
elseif any(diff(tspan(:)) <= 0) || size(tspan,2) < 2
  error('Field .tspan must be an increasing vector.');
end

% Check u0.
if ndims(umod.u0) == 2
  [Mspecies,Ncells] = size(umod.u0);
  Nreplicas = 1;
elseif ndims(umod.u0) == 3
  [Mspecies,Ncells,Nreplicas] = size(umod.u0);
else
  error('Incorrect number of dimensions of initial state.');
end
Ndofs = Mspecies*Ncells;
if issparse(umod.u0) || ~isa(umod.u0,'double')
  error('Initial state must be a double vector.');
elseif any(umod.u0(:) < 0)
  error('Initial state has negative elements.');
end

% Check D.
Ddiag = diag(umod.D);
if ~issparse(umod.D) || ~isa(umod.D,'double')
  error('Diffusion matrix must be sparse.');
elseif any(Ddiag > 0)
  error('Diffusion matrix has positive diagonal elements.');
elseif any(size(umod.D) ~= Ndofs)
  error('Wrong size of diffusion matrix.');
elseif any(abs(sum(umod.D,1)) > 1000*eps*abs(Ddiag)')
  error('Sum of rows criterion in diffusion matrix is violated.');
end

% Check N.
[MM,Mreactions] = size(umod.N);
if ~issparse(umod.N) || ~isa(umod.N,'double')
  error('Stochiometric matrix must be sparse.');
elseif MM ~= Mspecies
  error('Wrong size of stochiometric matrix.');
end

% Check G.
if ~issparse(umod.G) || ~isa(umod.G,'double')
  error('Dependency graph must be sparse.');
elseif any(size(umod.G) ~= [Mreactions,Mreactions+Mspecies])
  error('Wrong size of dependency graph.');
end

% Check vol.
vol = reshape(umod.vol,1,[]);
if issparse(vol) || ~isa(vol,'double')
  error('Volume must be a double vector.');
elseif any(vol <= 0)
  error('Volume vector must be positive.');
elseif size(vol,2) ~= Ncells
  error('Wrong size of volume vector.');
end

% Check sd.
sd = reshape(umod.sd,1,[]);
if issparse(sd) || ~isa(sd,'double')
  error('Subdomain must be a double vector.');
elseif size(sd,2) ~= Ncells
  error('Wrong size of subdomain vector.');
elseif any(sd ~= ceil(sd))
  error('Subdomain vector must be integer.');
end

% the following fields are added automatically by urdme if they are
% missing:

% Check report.
if ~isscalar(umod.report)
  error('Report must be a scalar.');
end

% Check seed.
seed = reshape(umod.seed,1,[]);
if issparse(seed) || ~isa(seed,'double')
  error('Seed must be a double vector.');
elseif size(seed,2) ~= 1 && size(seed,2) ~= Nreplicas
  error('Wrong size of seed vector.');
elseif any(seed ~= ceil(seed))
  error('Seed vector must be integer.');
end

% Check inline propensities.
if ~isstruct(umod.inline_propensities)
  error('Inline propensities must be specified as a struct.');
end
% fields assumed to exists:
K = umod.inline_propensities.K;
I = umod.inline_propensities.I;
S = umod.inline_propensities.S;
% (other fields will be rejected provided urdme parses umod)
M1 = 0;
if ~isempty(K) || ~isempty(I)
  M1 = size(K,2);
  if issparse(K) || ~isa(K,'double') || size(K,1) ~= 3 || M1 > Mreactions || ...
     issparse(I) || ~isa(I,'double') || any(size(I) ~= [3 M1])
    error('Format mismatch in inline propensities.');
  end
  if any(I(:) < 1 | Mspecies < I(:)) || any(I(:) ~= ceil(I(:)))
    error('Index out of bounds in inline propensity.');
  end
end
if ~isempty(S)
  if ~issparse(S) || ~isa(S,'double') || size(S,2) ~= M1
    error('Format mismatch in inline propensities.');
  end
  if any(S(:) ~= ceil(S(:)))
    error('Subdomain listed for inline propensities must be integers.');
  end
end

% Check ldata.
[ldsize,NN] = size(umod.ldata);
if issparse(umod.ldata) || ~isa(umod.ldata,'double')
  error('Local data matrix must be a double matrix.');
elseif NN ~= Ncells
  error('Wrong size of local data matrix.');
end

% Check gdata.
if issparse(umod.gdata) || ~isa(umod.gdata,'double')
  error('Global data matrix must be a double vector.');
end

% Check solverargs.
if ~isfield(umod,'solverargs');
  error('No field .solverargs.');
end
if ~iscell(umod.solverargs)
  error('Solver arguments must be a cell vector.');
end
umod.solverargs = reshape(umod.solverargs,1,[]);

% Check makeargs.
if ~isfield(umod,'makeargs');
  error('No field .makeargs.');
end
if ~iscell(umod.makeargs)
  error('Make arguments must be a cell vector.');
end
umod.makeargs = reshape(umod.makeargs,1,[]);