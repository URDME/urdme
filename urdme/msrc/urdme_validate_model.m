function urdme_validate_model(umod)
%URDME_VALIDATE_MODEL Validate URDME model.
%   URDME_VALIDATE_MODEL(UMOD) attempts to validate the model
%   expressed in UMOD and throws an error if a reaction which
%   decreases species counts is associated with a positive rate at one
%   of the domain boundaries.
%
%   The approach is to first compile the UDS-solver and then perform
%   trial evaluations of the propensities of the model. This is not
%   fool-proof but will catch most common mistakes.
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

% S. Engblom 2022-04-12

% parse and re-compile for the UDS-solver
umod = urdme(umod,'solver','uds', ...
             'parse',1,'compile',1,'solve',0);

% find all reactions which may decrease a molecular count
if ndims(umod.u0) == 2
  [Mspecies,Ncells] = size(umod.u0);
elseif ndims(umod.u0) == 3
  [Mspecies,Ncells,~] = size(umod.u0);
end
% loop: all reactions which decreases a count
ireactions = find(any(umod.N < 0,1));
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
      R = mexrhs(1,repmat(full(X),1,Ncells),size(umod.N,2), ...
                 umod.vol,umod.ldata,umod.gdata,umod.sd, ...
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