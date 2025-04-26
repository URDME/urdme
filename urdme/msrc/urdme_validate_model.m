function urdme_validate_model(umod,rd)
%URDME_VALIDATE_MODEL Validate URDME model.
%   URDME_VALIDATE_MODEL(UMOD) attempts to validate the model
%   expressed in UMOD and throws an error if the model seems to be
%   missformed.
%
%   URDME_VALIDATE_MODEL(UMOD,'r') investigates the reactions only and
%   throws an error if a reaction which decreases species counts is
%   associated with a positive rate at one of the domain boundaries.
%
%   The approach is to first compile the UDS-solver and then perform
%   trial evaluations of the propensities of the model. This is not
%   fool-proof but will catch most common mistakes.
%
%   URDME_VALIDATE_MODEL(UMOD,'d') investigates the diffusion-part
%   only and throws an error or displays a warning if the diffusion
%   matrix is incorrect.
%
%   URDME_VALIDATE_MODEL(UMOD,'rd') performs both of the above
%   checks. This is the default.
%
%   Example:
%     clear umod
%     umod = rparse([],{'X > k*X > @' ...
%                       'X+Y > k*X > @' ...
%                       '@ > k > Y'},{'X' 'Y'},{'k' 1},'test.c');
%     umod.vol = ones(1,5);
%     umod.tspan = [0 100];
%     umod.u0 = 10*ones(2,5);
%     umod.D = sparse(2*5,2*5);
%     umod.sd = ones(1,5);
%     urdme_validate_model(umod);
%
%   See also URDME_VALIDATE.

% S. Engblom 2024-04-11 (diffusion part added)
% S. Engblom 2022-04-12

if nargin == 1, rd = 'rd'; end

% initial check
umod = urdme(umod,'parse',1,'compile',0,'solve',0);
if ndims(umod.u0) == 2
  [Mspecies,Ncells] = size(umod.u0);
elseif ndims(umod.u0) == 3
  [Mspecies,Ncells,~] = size(umod.u0);
end

% reactions
if any(rd == 'r')
  % parse and re-compile for the UDS-solver
  umod = urdme(umod,'solver','uds', ...
               'modelname',['urdme_validate_model_' umod.modelname], ...
               'parse',1,'compile',1,'solve',0);
  mexrhs = str2func([umod.mexname '_mexrhs']);

  % find all reactions which may decrease a molecular count
  ireactions = find(any(umod.N < 0,1));
  % loop: all reactions which decreases a count
  for i = ireactions
    % loop: all species which are affected negatively
    jspecies = find(umod.N(:,i) < 0)';
    for j = jspecies
      % loop: all counts that could go negative
      kcount = -umod.N(j,i);
      for k = 0:kcount-1
        X = 10*ones(Mspecies,1); % (>> 0)
        % a single count k at species j which is affected negatively by
        % reaction i
        X(j) = k;
        R = mexrhs(umod.mexhash, ...
                   1,repmat(full(X),1,Ncells),size(umod.N,2), ...
                   umod.vol,umod.ldata,umod.gdata, ... 
                   umod.ldata_time,umod.gdata_time, ... 
                   umod.sd, ... 
                   umod.inline_propensities.K, ...
                   umod.inline_propensities.I, ...
                   umod.inline_propensities.S);
        R = reshape(R,[],Ncells);
        neg = find(any(R < 0,2),1,'first');
        if ~isempty(neg)
          error(sprintf('Rate #%d may be become negative.',neg));
        end
        R = R(i,:);
        if any(R ~= 0,2)
          error(sprintf('Reaction #%d may produce negative states.',i));
        end
      end
    end
  end
end

% diffusions
if any(rd == 'd')
  % check diffusion of individual species
  for i = 1:Mspecies
    D = umod.D(i:Mspecies:end,i:Mspecies:end);
    Ddiag = diag(D);
    % this should be caught by urdme_validate in the initial urdme-call:
    if any(Ddiag > 0)
      error(sprintf(['Diffusion matrix for species %d ' ...
                     'has positive diagonal elements.',i]));
    end
    if any(abs(sum(D,1)) > 1000*eps*abs(Ddiag)')
      % investigate: where did it diffuse?
      for j = [1:i-1 i+1:Mspecies]
        D_ = umod.D(j:Mspecies:end,i:Mspecies:end);
        if any(abs(sum(D_,1)) > 1000*eps)
          % this can be ok but you should know what you are doing:
          warning('urdme_validate_model:species_diffusion', ...
                  sprintf('Seems species #%d diffuses into species #%d.',i,j));
        end
      end
    end
  end
end
