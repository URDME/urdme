% URDME_VALIDATE
%
% Validates the fields in umod (sanity-checks). Since the core 
% simulators perform a very limited error-checking of the input, it is
% important to check data with this function before calling a simulation
% routine. 
%
% See also: urdme, rdme2mat.

function urdme_validate(umod)
 
% Check tspan.
if ~isfield(umod,'tspan');
    error('No field tspan in umod, check model file');
end

tspan = reshape(umod.tspan,1,[]);
if issparse(tspan)||~isa(tspan,'double') 
  error('Input time-span must be a double vector.');
elseif any(diff(tspan(:))<=0)||size(tspan,2)<2
  error('Input time-span must be an increasing vector.');
end

% Check u0.
if ~isfield(umod,'u0');
    error('No field u0 in umod, check model file');
end

[Mspecies,Ncells]=size(umod.u0);
Ndofs=Mspecies*Ncells;
if issparse(umod.u0)||~isa(umod.u0,'double')
  error('Initial state must be a double vector.');
elseif any(umod.u0(:)<0)
  error('Initial state has negative elements.');
end


% Check D.
if ~isfield(umod,'D');
    error('No field D in umod, check model file');
end

Ddiag=diag(umod.D);
if ~issparse(umod.D)||~isa(umod.D, 'double')
    error('Diffusion matrix must be sparse.');
elseif any(Ddiag>0)
    error('Diffusion matrix has positive diagonal elements.');
elseif any(size(umod.D)~=Ndofs)
    error('Wrong size of diffusion matrix.');
elseif any(abs(sum(umod.D,1))>1000*eps*abs(Ddiag)')
    error('Sum of rows criterion in diffusion matrix is violated.');
end

% Check N.
if ~isfield(umod,'N');
    error('No field N in umod, check model file');
end

[MM,Mreactions]=size(umod.N);
if ~issparse(umod.N)||~isa(umod.N, 'double')
  error('Stochiometric matrix must be sparse.');
elseif MM~=Mspecies
  error('Wrong size of stochiometric matrix.');
end

% Check G.
if ~isfield(umod,'G');
    error('No field G in umod, check model file');
end
if ~issparse(umod.G)||~isa(umod.G,'double')
  error('Dependency graph must be sparse.');
elseif any(size(umod.G)~=[Mreactions,Mreactions+Mspecies])
  error('Wrong size of dependency graph.');
end

% Check vol.
if ~isfield(umod,'vol');
    error('No field vol in umod, check model file');
end
vol=reshape(umod.vol,1,[]);
if issparse(vol)||~isa(vol,'double')
  error('Volume must be a double vector.');
elseif any(vol<=0)
  error('Volume vector must be positive.');
elseif size(vol,2)~=Ncells
  error('Wrong size of volume vector.');
end

% Check sd.
if ~isfield(umod,'sd');
    error('No field sd in umod, check model file');
end
sd=reshape(umod.sd,1,[]);
if issparse(sd)||~isa(sd,'double')
  error('Subdomain must be a double vector.');
elseif size(sd,2)~=Ncells
  error('Wrong size of subdomain vector.');
elseif any(sd~=ceil(sd))
  error('Subdomain vector must be integer.');
end

% Check data.
if ~isfield(umod,'data');
    error('No field data in umod, check model file');
end
[dsize,NN]=size(umod.data);
if issparse(umod.data)||~isa(umod.data,'double')
   error('Data matrix must be a double matrix.');
elseif NN~=Ncells
   error('Wrong size of data matrix.');
end

% check OPTS.
% if( size(OPTS,2) < 4)
%   error(['Option vector must be at least 4 elements. Are you defining ' ...
%          'all default parameters?']);
% elseif issparse(OPTS)||~isa(OPTS,'double')
%   error('Options must be a double vector.');
% end

% % Check K.
% if issparse(K)||~isa(K,'double')
%   error('Inline propensity coefficient must be a double matrix.');
% elseif size(K,1)~=3
%   error('Wrong size of inline propensity coefficient matrix.');
% end
% 
% % Check I.
% if issparse(I)||~isa(I,'double')
%   error('Inline propensity index must be a double matrix.');
% elseif any(size(I)~=size(K)) || size(I,2)>Mreactions
%   error('Wrong size of inline propensity index matrix.');
% elseif any(I(:)<1 | Mspecies<I(:))
%   error('Species index out of bounds.');
% end
% 
% % Check S.
% if ~issparse(S)||~isa(S,'double')
%   error('Inline subdomain matrix must be sparse.');
% elseif size(S,2)~=size(I,2)
%   error('Wrong size of inline subdomain matrix.');
% end

