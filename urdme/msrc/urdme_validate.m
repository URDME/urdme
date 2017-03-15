function urdme_validate(umod)
% URDME_VALIDATE Validate URDME struct.
%    URDME_VALIDATE(UMOD) attempts to validate the fields in UMOD and
%    throws an error if a field is missing or has an incorrect format.
%
%    See also URDME.

% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% J. Cullhed 2008-06-18 (rdme.m)

% Check tspan.
if ~isfield(umod,'tspan');
  error('No field .tspan.');
end
tspan = reshape(umod.tspan,1,[]);
if issparse(tspan) || ~isa(tspan,'double')
  error('Field .tspan must be a double vector.');
elseif any(diff(tspan(:)) <= 0) || size(tspan,2) < 2
  error('Field .tspan must be an increasing vector.');
end

% Check u0.
if ~isfield(umod,'u0');
  error('No field .u0.');
end
[Mspecies,Ncells] = size(umod.u0);
Ndofs = Mspecies*Ncells;
if issparse(umod.u0) || ~isa(umod.u0,'double')
  error('Initial state must be a double vector.');
elseif any(umod.u0(:) < 0)
  error('Initial state has negative elements.');
end

% Check D.
if ~isfield(umod,'D');
  error('No field .D.');
end

Ddiag = diag(umod.D);
if ~issparse(umod.D) || ~isa(umod.D,'double')
  error('Diffusion matrix must be sparse.');
elseif any(Ddiag > 0)
  error('Diffusion matrix has positive diagonal elements.');
elseif any(size(umod.D) ~= Ndofs)
  error('Wrong size of diffusion matrix.');
elseif any(abs(sum(umod.D,1))>1000*eps*abs(Ddiag)')
  error('Sum of rows criterion in diffusion matrix is violated.');
end

% Check N.
if ~isfield(umod,'N');
  error('No field .N.');
end

[MM,Mreactions] = size(umod.N);
if ~issparse(umod.N) || ~isa(umod.N,'double')
  error('Stochiometric matrix must be sparse.');
elseif MM ~= Mspecies
  error('Wrong size of stochiometric matrix.');
end

% Check G.
if ~isfield(umod,'G');
    error('No field .G.');
end
if ~issparse(umod.G) || ~isa(umod.G,'double')
  error('Dependency graph must be sparse.');
elseif any(size(umod.G) ~= [Mreactions,Mreactions+Mspecies])
  error('Wrong size of dependency graph.');
end

% Check vol.
if ~isfield(umod,'vol');
  error('No field .vol.');
end
vol = reshape(umod.vol,1,[]);
if issparse(vol) || ~isa(vol,'double')
  error('Volume must be a double vector.');
elseif any(vol <= 0)
  error('Volume vector must be positive.');
elseif size(vol,2) ~= Ncells
  error('Wrong size of volume vector.');
end

% Check sd.
if ~isfield(umod,'sd');
  error('No field .sd.');
end
sd = reshape(umod.sd,1,[]);
if issparse(sd) || ~isa(sd,'double')
  error('Subdomain must be a double vector.');
elseif size(sd,2) ~= Ncells
  error('Wrong size of subdomain vector.');
elseif any(sd ~= ceil(sd))
  error('Subdomain vector must be integer.');
end

% the following fields are added automatically bu urdme if they are
% missing:

% Check ldata.
if ~isfield(umod,'ldata');
  error('No field .ldata.');
end
[ldsize,NN] = size(umod.ldata);
if issparse(umod.ldata) || ~isa(umod.ldata,'double')
  error('Local data matrix must be a double matrix.');
elseif NN ~= Ncells
  error('Wrong size of local data matrix.');
end

% Check gdata.
if ~isfield(umod,'gdata');
  error('No field .gdata.');
end
if issparse(umod.gdata) || ~isa(umod.gdata,'double')
  error('Global data matrix must be a double vector.');
end

% Check solverargs.
if ~isfield(umod,'solverargs');
  error('No field .solverargs.');
end
