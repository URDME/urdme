function [K,I,N,G] = rparse_inline(umod,r,spec,rate,seq)
%RPARSE_INLINE URDME inline reaction propensity parser.
%   UMOD = RPARSE_INLINE(UMOD,R,SPEC,RATE,SEQ) creates/augments the
%   URDME structure UMOD with fields describing the reactions R with
%   species SPEC, and rate constants RATE. The sequence expansion
%   argument SEQ is optional, see functionality below.
%
%   Specific syntaxes and descriptions now follow. Note: adding a
%   (possible empty) URDME structure UMOD to the inputs as in the main
%   call above implies that this structure is
%   created/augmented. Otherwise the results are returned in the
%   various outputs as specified below.
%
%   [K,I] = RPARSE_INLINE(R,SPEC,RATE) constructs inline propensities
%   (K,I) for the reactions R with species SPEC and rate-constants
%   RATE.
%
%   R is a cell-vector containing reactions on the form 'X+Y+... >
%   rate constant k > U+V+...' The left (right) side is the initial
%   (final) state and the rate constant is written in between the
%   '>'-signs. The special symbol '@' is reserved for the empty
%   set. For example,
%     R = {'X+Y > mu > Z' 'Z > k > X+Y'}
%   expresses a reversible bimolecular reaction. Note that, for inline
%   propensities, the combinatorial law and the scaling with volume is
%   implicitly understood. Hence with RPARSE, the above is equivalent
%   to
%     R = {'X+Y > mu*X*Y/vol > Z' 'Z > k*Z > X+Y'}.
%   Note also the implicit scaling with 2 for the case of a
%   dimerization: R = {'X+X > nu > @'} with RPARSE_INLINE is
%   equivalent to R = {'X+X > nu*X*(X-1)/(2*vol) > @'} with RPARSE,
%   see NSM for more details.
%
%   SPEC contains the names of the involved species, for example,
%     SPEC = {'X' 'Y' 'Z'}.
%
%   RATE contains the names of the rate-constants followed by their
%   values, for example,
%     RATE = {'k' 1 'mu' 0.5e-7}.
%
%   [K,I,N,G] = RPARSE_INLINE(...) additionally outputs the
%   stoichiometric matrix N and the dependency graph G, to be placed
%   in the fields umod.N and umod.G of the URDME structure.
%
%   Examples:
%     % linear birth-death example
%     [K1,I1] = rparse_inline({'@ > k > X' 'X > mu > @' }, ...
%                             {'X'},{'k' 1 'mu' 1e-3});
%
%     % 2 reacting species
%     [K2,I2,N,G] = rparse_inline({'@ > k > X' ...
%                                  '@ > k > Y' ...
%                                  'X > mu > @' ...
%                                  'Y > mu > @' ...
%                                  'X+Y > kk > @'}, ...
%                                  {'X' 'Y'},{'k' 1 'mu' 1e-3 'kk' 1e-4});
%
%   See also NSM, RPARSE, URDME.

% S. Engblom 2019-11-26 (Revision, now creating/augmenting URDME structure)
% S. Engblom 2017-03-01

% umod-in-umod-out syntax?
if ~iscell(umod) && nargout <= 1
  % rparse(umod,r,spec,rate,[seq])
  if ~isstruct(umod)
    if ~isempty(umod)
      error('URDME structure expected as first argument.');
    end
    umod = struct; % (rparse used as constructor)
  end
  if nargin < 4
    error('Minimum 4 arguments required to create/augment URDME structure.');
  end
else
  % rparse(umod,r,spec,rate,seq) --> rparse(r,spec,rate,[seq])
  if nargin > 3
    seq = rate;
  end
  rate = spec;
  spec = r;
  r = umod;
  umod = [];
end
% from now on, just ignore umod until the end

r = reshape(r,1,[]);
spec = reshape(spec,1,[]);
rate = reshape(rate,1,[]);
ratedef = rate(2:2:end);
rate = rate(1:2:end);

% some checks
if ~all(cellfun('isclass',spec,'char'))
  error('Species must be specified in a cell-vector of strings.');
elseif any(size(rate) ~= size(ratedef)) || ...
      ~all(cellfun('isclass',rate,'char')) || ...
      ~all(cellfun('isclass',ratedef,'double'))
  error('Rates must be specified as property/value-pairs.');
elseif size(unique(spec),2) ~= size(spec,2) || ...
      size(unique(rate),2) ~= size(rate,2)
  error('Species and rates must consist of unique names.');
elseif ~isempty(intersect(spec,rate))
  error('Species and rates must use different names.');
end

% sequence expansion
if nargin > (3+isstruct(umod))
  r = seqexpand(r,seq);
  spec = seqexpand(spec,seq);
  [rate,ratedef] = seqexpand(rate,ratedef,seq);
end

K = zeros(3,size(r,2));
I = ones(3,size(r,2));
N = sparse(size(spec,2),size(r,2));
H = sparse(size(spec,2),size(r,2));
for i = 1:size(r,2)
  % separate reaction r{i} according to 'from > propensity > dest'
  ixgt1 = find(r{i} == '>',1);
  ixgt2 = find(r{i} == '>',1,'last');
  if isempty(ixgt1) || isempty(ixgt2) || ixgt1 == ixgt2
    error('Each reaction must contain at least 2 ''>''.');
  end
  from = r{i}(1:ixgt1-1);
  prop = r{i}(ixgt1+1:ixgt2-1);
  dest = r{i}(ixgt2+1:end);

  % remove spaces and the empty set
  from = from(from ~= ' ' & from ~= '@');
  prop = prop(prop ~= ' ');
  dest = dest(dest ~= ' ' & dest ~= '@');

  % split from and dest into 'spec1 + spec2 + ..'
  from = l_split(from,'+');
  dest = l_split(dest,'+');

  % assign each species its number according to the ordering in spec
  [foo,ifrom] = ismember(from,spec);
  [foo,idest] = ismember(dest,spec);
  if any(ifrom == 0)
    error(['Unknown species ''' char(from{find(ifrom == 0,1,'first')}) ...
           ''' in reaction #' num2str(i) '.']);
  end
  if numel(ifrom) > 2
    error(['Reaction #' num2str(i) ' is not elementary.']);
  end
  if any(idest == 0)
    error(['Unknown species ''' char(dest{find(idest == 0,1,'first')}) ...
           ''' in reaction #' num2str(i) '.']);
  end

  % the corresponding column of the stoichiometric matrix is now known
  N(:,i) = sparse(ifrom,1,-1,size(spec,2),1)+ ...
           sparse(idest,1,1,size(spec,2),1);

  % the contribution to the dependency graph is now known
  H(:,i) = sparse(ifrom,1,1,size(spec,2),1);

  % determine the rate constant
  [foo,irate] = ismember(prop,rate);
  if irate == 0
    error(['Unknown rate in reaction #' num2str(i) '.']);
  end
  K(3-numel(ifrom),i) = ratedef{irate};

  % index of involved species
  switch numel(ifrom)
   case 2, I(1:2,i) = ifrom;
   case 1, I(3,i) = ifrom;
  end
end

% after finding out the stoichiometric matrix we can converge on G
G = [H' H'*abs(N)]; % [diffusion reaction]-parts of G
G = double(G ~= 0);

% umod-in-umod-out syntax?
if isstruct(umod)
  % output: UMOD.{inline_propensities.{K I S} N G 
  % private.{Reactions Species RateNames RateVals}}
  umod.inline_propensities.K = K;
  umod.inline_propensities.I = I;
  umod.inline_propensities.S = sparse(0,size(K,2));
  umod.N = N;
  umod.G = G;
  umpd.private.Reactions = r;
  umod.private.Species = spec;
  umod.private.RateNames = rate;
  umod.private.RateVals = ratedef;
  K = umod;
end

%---------------------------------------------------------------------------
function c = l_split(s,d)
%L_SPLIT Split string.
%   C = L_SPLIT(S,D) splits the string S at each occurrence of the character
%   D and returns all intermediate strings in the cell-vector C.
%  
%   Example:
%     l_split('A+B/C + D','+') returns {'A' 'B/C ' ' D'}

% handle empty case
if isempty(s), c = {}; return; end

ii = [0 find(s == d) size(s,2)+1];
c = cell(1,size(ii,2)-1);
for i = 1:size(ii,2)-1
  c{i} = s(ii(i)+1:ii(i+1)-1);
end

%---------------------------------------------------------------------------
