function umod = rparse(umod,r,spec,rate,filename,seq)
%RPARSE URDME reaction propensity parser.
%   UMOD = RPARSE(UMOD,R,SPEC,RATE,FILENAME,SEQ) creates/augments the
%   URDME structure UMOD with fields describing the reactions R with
%   species SPEC, rate constants RATE, and propensity file
%   FILENAME. If the input UMOD is empty the structure is
%   simultaneously created. The sequence expansion argument SEQ is
%   optional, see functionality below.
%
%   R is a cell-vector containing reactions on the form
%      'X+Y+... > F(...) > U+V+...'
%   The left (right) side is the initial (final) state and the propensity is
%   written in between the '>'-signs. The special symbol '@' is reserved for
%   the empty set. For example,
%     R = {'X+Y > mu*X*Y > Z' 'Z > k*Z > X+Y'}
%   expresses a reversible bimolecular reaction.
%
%   SPEC contains the names of the involved species, for example,
%     SPEC = {'X' 'Y' 'Z'}.
%
%   RATE contains the names of the rate constants followed by their
%   values, for example,
%     RATE = {'k' 1 'mu' 0.5e-7}.
%   Rate constants can be vectors, for example,
%     RATE = {'birth' [1 0.5 0.2]},
%   which implies that the rate constants birth[0..2] can be used in
%   the propensity expressions. Furthermore, you can use 'gdata' and
%   'ldata' to denote rate constants to be passed via these fields,
%   for example,
%     RATE = {'alpha' 'gdata' 'beta' 'ldata'},
%   which will reauire that the fields umod.gdata and umod.data are
%   initialized accordingly. A yet more advanced use is the
%   time-dependent fields, for example,
%     RATE = {'alpha' 'gdata_time' 'beta' 'gdata_time'}.
%   With umod.data_time = [0 5 10] and umod.gdata_time a 2-by-3
%   matrix, this changes the parameters alpha and beta at the given
%   points in time.
%
%   FILENAME is the name of the file to write the result to; a string,
%   stored in the field umod.propensity. If no extension is given,
%   '.c' is added. This file, if already existing, is only overwritten
%   provided that is was originally created by RPARSE. The resulting
%   C-file is ready to be compiled and linked with any of the solvers
%   in URDME. If FILENAME is empty or missing, the generated code is
%   written out explicitly (set FILENAME = '~' to skip this).
%
%   SEQ allows for sequence expansion as determined by SEQ. See
%   SEQEXPAND for further details of this powerful syntax.
%
%   The field UMOD.private.rp contains the various inputs to RPARSE as
%   well as some additionally created information.
%
%   UMOD = RPARSE(UMOD,'jac') is a specialized syntax for
%   Jacobians. It is assumed that UMOD has been obtained via an
%   initial call to RPARSE. The second call then augments the
%   propensity file with a code for the Jacobian, see MEXJAC.
%
%   Examples:
%     % propensity code for linear birth-death model
%     umod = rparse([],{'@ > k*vol > X' 'X > mu*X > @' }, ...
%                      {'X'},{'k' 1 'mu' 1e-3});
%
%     % 2 reacting species, create URDME structure
%     vmod = rparse([],{'@ > k*vol > X' ...
%                       '@ > k*vol > Y' ...
%                       'X > mu*X > @' ...
%                       'Y > mu*Y > @' ...
%                       'X+Y > kk*X*Y/vol > @'}, ...
%                       {'X' 'Y'},{'k' 1 'mu' 1e-3 'kk' 1e-4}, ...
%                       'bimol');
%
%   See also RPARSE_INLINE, SEQEXPAND, URDME, NSM.

% S. Engblom 2024-06-14 (removed old calling sequence)
% S. Engblom 2024-06-10 ('jac'-syntax)
% S. Engblom 2024-05-08 (Revision, ldata_time/gdata_time)
% S. Engblom 2024-05-02 (Revision, ghost species/rates)
% S. Engblom 2019-11-25 (Revision, creating/augmenting umod, ldata/gdata-syntax)
% S. Engblom 2019-11-06 (Revision, now using URDMEstate_t)
% S. Engblom 2017-0-08 (input seq added)
% S. Engblom 2017-03-01 (URDME 1.3)
% S. Engblom 2016-01-07 (empty default rate)
% S. Engblom 2007-04-06

% syntax
if nargout > 1
  error('One output only.');
end
if nargin == 2 && strcmp(r,'jac')
  umod = l_rparse_jac(umod);
  return;
end
if nargin < 6
  seq = [];
  if nargin < 5
    filename = '';
    if nargin < 4
      error('Minimum 4 input arguments required.');
    end
    
  end
end
if ~isstruct(umod)
  if ~isempty(umod)
    error('URDME structure expected as first argument.');
  end
  umod = struct; % (rparse used as constructor)
end

r = reshape(r,1,[]);
spec = reshape(spec,1,[]);
% also allow rates defined by a struct
if isstruct(rate)
  ratedef = struct2cell(rate)';
  rate = fieldnames(rate)';
else
  rate = reshape(rate,1,[]);
  ratedef = rate(2:2:end);
  rate = rate(1:2:end);
end

% some checks
if ~all(cellfun('isclass',spec,'char'))
  error('Species must be specified in a cell-vector of strings.');
elseif any(size(rate) ~= size(ratedef)) || ...
      ~all(cellfun('isclass',rate,'char')) || ...
      ~all(cellfun('isclass',ratedef,'double') | ...
           cellfun('isclass',ratedef,'char'))
  error('Rates must be specified as property/value-pairs.');
elseif size(unique(spec),2) ~= size(spec,2) || ...
      size(unique(rate),2) ~= size(rate,2)
  error('Species and rates must consist of unique names.');
elseif ~isempty(intersect(spec,rate))
  error('Species and rates must use different names.');
end
reserved = {'xstate' 'time' 'vol' 'ldata' 'gdata' ...
            'ldata_time' 'gdata_time' 'sd'};
res1 = intersect(spec,reserved);
res2 = intersect(rate,reserved);
if ~isempty(res1)
  error(sprintf(['Species should not clash with propensity ' ...
                 'input arguments: ''%s''.'],res1{1}));
elseif ~isempty(res2)
  error(sprintf(['Rates should not clash with propensity ' ...
                 'input arguments: ''%s''.'],res2{1}));
end

% sequence expansion
if nargin > (4+isstruct(umod))
  r = seqexpand(r,seq);
  spec = seqexpand(spec,seq);
  [rate,ratedef] = seqexpand(rate,ratedef,seq);
end

N = sparse(size(spec,2),size(r,2));
H = sparse(size(spec,2),size(r,2));
D = sparse(size(rate,2),size(r,2));
prop = cell(size(r));
for i = 1:size(r,2)
  % separate reaction r{i} according to 'from > propensity > dest'
  ixgt1 = find(r{i} == '>',1);
  ixgt2 = find(r{i} == '>',1,'last');
  if isempty(ixgt1) || isempty(ixgt2) || ixgt1 == ixgt2
    error('Each reaction must contain at least 2 ''>''.');
  end
  from = r{i}(1:ixgt1-1);
  prop{i} = r{i}(ixgt1+1:ixgt2-1);
  dest = r{i}(ixgt2+1:end);

  % remove spaces and the empty set
  from = from(from ~= ' ' & from ~= '@');
  prop{i} = prop{i}(prop{i} ~= ' ');
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
  if any(idest == 0)
    error(['Unknown species ''' char(dest{find(idest == 0,1,'first')}) ...
           ''' in reaction #' num2str(i) '.']);
  end

  % the corresponding column of the stoichiometric matrix is now known
  N(:,i) = sparse(ifrom,1,-1,size(spec,2),1)+ ...
           sparse(idest,1,1,size(spec,2),1);

  % search for species and rates and rewrite in intermediate format
  [rprop{i},deps,depr] = l_rewriteprop(prop{i},spec,rate);

  % the contribution to the dependency graph is now known
  H(:,i) = sparse(deps,1,1,size(spec,2),1);

  % accumulate also dependecies on rates
  D(:,i) = sparse(depr,1,1,size(rate,2),1);
end

% after finding out the stoichiometric matrix we can converge on G
G = [H' H'*abs(N)]; % [diffusion reaction]-parts of G

% reactions dependent on time-dependent rates
tdrates = D(strcmp(ratedef,'ldata_time') | strcmp(ratedef,'gdata_time'),:);
if ~isempty(tdrates)
  G = [G any(tdrates,1)']; % (one more column to the right)
end
G = double(G ~= 0);

% suspicious:
ghost_species = find(~(any(H,2) | any(N,2)));
if ~isempty(ghost_species)
  ghosts = sprintf('{ %s}',sprintf('#%d ',ghost_species));
  warning('rparse:ghost_species', ...
          ['One or more species ' ghosts ' do not participate in any reactions ' ...
           'and are not affected by any reactions.']);
end
ghost_rates = find(~any(D,2));
if ~isempty(ghost_rates)
  ghosts = sprintf('{ %s}',sprintf('#%d ',ghost_rates));
  warning('rparse:ghost_rates', ...
          ['One or more rates ' ghosts ' do not participate in any reactions.']);
end

% heading
F = ['/* [Remove/modify this line not to overwrite this file] */\n'];
F = [F '/* Generated by RPARSE ' datestr(now,'yyyy-mm-dd HH:MM') ' */\n\n'];
F = [F '/* Reactions:\n'];
for i = 1:size(r,2)
  F = [F '     ' char(r{i}) '\n'];
end
F = [F '*/\n\n'];

% includes
F = [F '#include <math.h>\n' ...
     '#include "propensities.h"\n' '#include "report.h"\n\n'];

% species enum
if size(spec,2) > 0
  F = [F 'enum Species {'];
  for i = 1:size(spec,2)
    F = [F sprintf('\n  %s,',char(spec{i}))];
  end
  F(end) = [];
  F = [F '\n};\n\n'];
end

% number of reactions
F = [F sprintf('const int NR = %d; /* number of reactions */\n\n', ...
               size(N,2))];

% rate constants
F = [F '/* rate constants */\n'];
ldenum = ''; ldrate = {};
gdenum = ''; gdrate = {};
ldtenum = ''; ldtrate = {};
gdtenum = ''; gdtrate = {};
Frate = '';
for i = 1:size(rate,2)
  % local/global/time-dependent data?
  if ischar(ratedef{i})
    if strcmp(ratedef{i}(:)','ldata')
      % (note: seqexpand produces column strings)
      if isempty(ldenum)
        ldenum = 'enum ldataRATE {';
      end
      ldenum = [ldenum sprintf('\n  %s,',char(rate{i}))];
      ldrate = [ldrate rate{i}];
    elseif strcmp(ratedef{i}(:)','gdata')
      if isempty(gdenum)
        gdenum = 'enum gdataRATE {';
      end
      gdenum = [gdenum sprintf('\n  %s,',char(rate{i}))];
      gdrate = [gdrate rate{i}];
    elseif strcmp(ratedef{i}(:)','ldata_time')
      if isempty(ldtenum)
        ldtenum = 'enum ldata_timeRATE {';
      end
      ldtenum = [ldtenum sprintf('\n  %s,',char(rate{i}))];
      ldtrate = [ldtrate rate{i}];
    elseif strcmp(ratedef{i}(:)','gdata_time')
      if isempty(gdtenum)
        gdtenum = 'enum gdata_timeRATE {';
      end
      gdtenum = [gdtenum sprintf('\n  %s,',char(rate{i}))];
      gdtrate = [gdtrate rate{i}];
    else
      error(['String rates must be one of {''ldata'',''gdata'',' ...
             '''ldata_time'',''gdata_time''}.']);
    end
  % the magical constant 21 is supposed to be the DECIMAL_DIG of float.h
  elseif isscalar(ratedef{i})
    Frate = sprintf('%.*e',21,ratedef{i});
    % safety: read the constant back and see that it is the same
    if abs(sscanf(Frate,'%e')-ratedef{i}) > eps(1000)*abs(ratedef{i})
      warning(sprintf(['Large numerical errors in character conversion ' ...
                       'of rate. Conversion result: %s.'],Frate));
    end
    F = [F sprintf('const double %s = ',rate{i}) Frate ';\n'];
  else
    % note: the user is responsible for correctly indexing the elements of
    % the rate-vector
    F = [F sprintf('const double %s[] = {',rate{i})];
    for j = 1:numel(ratedef{i})
      Frate = sprintf('%.*e',21,ratedef{i}(j));
      if abs(sscanf(Frate,'%e')-ratedef{i}(j)) > eps(1000)*abs(ratedef{i}(j))
        warning(sprintf(['Large numerical errors in character conversion ' ...
                         'of rate. Conversion result: %s.'],Frate));
      end
      F = [F Frate ','];
    end
    F = [F(1:end-1) '};\n'];
  end
end
if ~isempty(Frate)
  F = [F '\n'];
end
if ~isempty(ldenum)
  ldenum(end) = [];
  ldenum = [ldenum '\n};\n\n'];
  F = [F ldenum];
end
if ~isempty(gdenum)
  gdenum(end) = [];
  gdenum = [gdenum '\n};\n\n'];
  F = [F gdenum];
end
if ~isempty(ldtenum)
  ldtenum(end) = [];
  ldtenum = [ldtenum '\n};\n\n'];
  F = [F ldtenum];
end
if ~isempty(gdtenum)
  gdtenum(end) = [];
  gdtenum = [gdtenum '\n};\n\n'];
  F = [F gdtenum];
end

% finish propensities
for i = 1:size(r,2)
  % search for species and rates in the propensity, rewrite as
  % 'xstate[<species>]', 'ldata[<rate>]', 'gdata[<rate>]',
  % 'ldata_time[<rate>]', 'gdata_time[<rate>]', and keep the rest of
  % the rates as constants
  rprop{i} = l_finishprop(rprop{i},spec,rate,ldrate,gdrate,ldtrate,gdtrate);
end

% declaration of propensity functions
F = [F '/* forward declarations */\n'];
for i = 1:size(r,2)
  F = [F sprintf(['double rFun%d(const URDMEstate_t *xstate,' ...
                  'double time,double vol,\n' ...
                  '             const double *ldata,' ...
                  'const double *gdata,\n' ...
                  '             const double *ldata_time,' ...
                  'const double *gdata_time,int sd);\n'],i)];
end
F = [F '\n'];

% static propensity vector
F = [F '/* static propensity vector */\n'];
F = [F 'static PropensityFun ptr[] = {'];
for i = 1:size(r,2)
  F = [F sprintf('rFun%d,',i)];
end
F(end) = [];
F = [F '};\n\n'];

% definition of propensities
F = [F '/* propensity definitions */\n'];
for i = 1:size(r,2)
  F = [F sprintf(['double rFun%d(const URDMEstate_t *xstate,' ...
                  'double time,double vol,\n' ...
                  '             const double *ldata,' ...
                  'const double *gdata,\n' ...
                  '             const double *ldata_time,' ...
                  'const double *gdata_time,int sd)\n{\n'],i)];
  F = [F sprintf('  return %s;\n}\n\n',rprop{i})];
end

% allocation/deallocation
F = [F '/* URDME solver interface */\n'];
F = [F 'PropensityFun *ALLOC_propensities(size_t Mreactions)\n'];
F = [F '{\n'];
F = [F '  if (Mreactions > NR) PERROR("Wrong number of reactions.");\n'];
F = [F '  return ptr;\n}\n\n'];

F = [F 'void FREE_propensities(PropensityFun *ptr)\n'];
F = [F '{ /* do nothing since a static array was used */ }\n\n'];

% write to file
if ~isempty(filename) && ~isscalar(filename) && all(filename ~= '~')
  % add extension if missing
  if ~any(filename == '.'), filename = [filename '.c']; end

  % check if file already exists
  if exist(filename,'file')
    [fid,msg] = fopen(filename,'rt');
    if fid == -1, error(msg); end
    s = fgets(fid);
    if ~strncmp(s,'/* [Remove/modify this line not to overwrite this file] */',58)
      warning('Will not overwrite existing file if not created by RPARSE.');
      filename = '';
    end
    fclose(fid);
  end

  % write to file
  if ~isempty(filename)
    [fid,msg] = fopen(filename,'wt');
    if fid == -1, error(msg); end
    fprintf(fid,F);
    fclose(fid);
  end
elseif any(filename ~= '~')
  % dump code to screen
  fprintf(1,F);
end
F = sprintf(F);

% LaTeX output
if nargout > 3 || isstruct(umod)
  L = ['\\begin{align}\n' ...
       '  \\left. \\begin{array}{rcl}\n'];
  for i = 1:size(r,2)
    ixgt = find(r{i} == '>');
    from_ = r{i}(1:ixgt(1)-1);
    prop_ = r{i}(ixgt(1)+1:ixgt(2)-1);
    dest_ = r{i}(ixgt(2)+1:end);
    from_ = strrep(from_,'@','\\emptyset');
    dest_ = strrep(dest_,'@','\\emptyset');
    L = [L '    ' from_ ' & \\xrightarrow{' prop_ '} & ' dest_ '\t\\\\\n'];
  end

  L = [L '  \\end{array} \\right\\}.\n' ...
       '\\end{align}'];
  L = sprintf(L);
end

% output: UMOD.{propensities N G <ldata> <gdata> 
% private.rp.{Reactions Species RateNames RateVals Propensities LaTeX}}
umod.propensities = filename;
% augment non-empty structure?
if isfield(umod,'N') || isfield(umod,'G')
  % when augmenting; by convention, put rparse's propensities last:
  try
    Nr1 = size(umod.N,2);
    Nr2 = size(N,2);
    umod.N = [umod.N N];
    Nspecies = size(umod.N,1);
    % how old/new (Nr1/Nr2) reactions affect new/old (Nr2/Nr1) ones we
    % cannot know without parsin, so we have to assume a full
    % coupling:
    sp1 = sparse(ones(Nr2,Nr1));
    umod.G = [[umod.G(:,1:Nspecies); G(:,1:Nspecies)] ...
              [[umod.G(:, Nspecies+1:end); sp1] ...
               [sp1'; G(:,Nspecies+1:end)]]];
  catch
    error(['Unable to augment URDME structure; existing species ' ...
           'seem to mismatch.']);
  end
else
  umod.N = N;
  umod.G = G;
end
if ~isempty(ldrate), umod.ldata = ldrate; end
if ~isempty(gdrate), umod.gdata = gdrate; end
if ~isempty(ldtrate), umod.ldata_time = ldtrate; end
if ~isempty(gdtrate), umod.gdata_time = gdtrate; end
% save rparse-description in "rp"-struct:
umod.private.rp.Reactions = r;
umod.private.rp.Species = spec;
umod.private.rp.RateNames = rate;
umod.private.rp.RateVals = ratedef;
% parsed/partially parsed items:
umod.private.rp.prop = prop;
umod.private.rp.rprop = rprop;
umod.private.rp.ldrate = ldrate;
umod.private.rp.gdrate = gdrate;
umod.private.rp.ldtrate = ldtrate;
umod.private.rp.gdtrate = gdtrate;
% "bonus":
umod.private.rp.LaTeX = L;
umod.private.rp.code = F;

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
function [prop,deps,depr] = l_rewriteprop(prop,spec,rate)
%L_REWRITEPROP Rewrite propensity.
%   [RPROP,DEPS,DEPR] = L_REWRITEPROP(PROP,SPEC,RATE) rewrites the
%   propensity PROP by replacing all species by '$[j]' where j is the
%   numbering in SPEC, and all rates by '$[-j]$' where j is the
%   numbering in RATE. On return, DEPS/DEPR contain all species/rates
%   upon which the propensity depends.

% The intermediate representation is '$[<name>]'. We replace larger
% strings first in order to avoid exchanging e.g., 'BA' with '$[1]A'
% whenever 'B' and 'BA' are two different species/rates
symbols = [spec; ...
           mat2cell(1:size(spec,2),1,ones(1,size(spec,2)))];
% rates have a negative sign:
symbols = [symbols [rate; ...
                    mat2cell(-(1:size(rate,2)),1,ones(1,size(rate,2)))]];
[foo,i] = sort(cellfun('prodofsize',symbols(1,:)));
symbols = symbols(:,i);
deps = zeros(1,0);
depr = zeros(1,0);
for j = size(symbols,2):-1:1
  if strfind(prop,symbols{1,j})
    if symbols{2,j} > 0
      % dependent on species
      deps = [deps symbols{2,j}];
    else
      % dependent on rate
      depr = [depr -symbols{2,j}];
    end
    % replace with $[<name>]
    prop = strrep(prop,symbols{1,j},['$[' num2str(symbols{2,j}) ']']);
  end
end

%---------------------------------------------------------------------------
function prop = l_finishprop(prop,spec,rate,ldrate,gdrate,ldtrate,gdtrate)
%L_FINISHPROP Finish propensities.
%   PROP = L_FINISHPROP(PROP,SPEC,RATE,LDRATE,GDRATE,LDTRATE,GDTRATE))
%   replaces back the species/rates names from the intermediate
%   format. In the C-code, enums are used for the species, and for the
%   ldata/gdata/ldata_time/gdata_time rate constants, while for the
%   other rate constants this is just a declared constant.
offset = 1;
while 1
  % find prop(i:j) = '$[<name>]'
  i = find(prop(offset:end) == '$',1);
  if isempty(i), break; end
  i = i+offset-1;
  j = find(prop(i+2:end) == ']',1);
  j = j+i+2-1;
  n = str2num(prop(i+2:j-1));
  if n > 0
    % species $ --> $1[...]
    prop = [prop(1:i) '1[' spec{n} prop(j:end)];
    offset = i+3;
  else
    % rates: ldata --> $2[...], gdata --> $3[...], ldata_time --> $4[...],
    % gdata_time --> $5[...], otherwise just remove the '$[]'
    if any(strcmp(rate{-n},ldrate))
      prop = [prop(1:i) '2[' rate{-n} prop(j:end)];
      offset = i+3;
    elseif any(strcmp(rate{-n},gdrate))
      prop = [prop(1:i) '3[' rate{-n} prop(j:end)];
      offset = i+3;
    elseif any(strcmp(rate{-n},ldtrate))
      prop = [prop(1:i) '4[' rate{-n} prop(j:end)];
      offset = i+3;
    elseif any(strcmp(rate{-n},gdtrate))
      prop = [prop(1:i) '5[' rate{-n} prop(j:end)];
      offset = i+3;
    else
      prop = [prop(1:i-1) rate{-n} prop(j+1:end)];
      offset = i;
    end
  end
end

% final replace
prop = strrep(prop,'$1','xstate');
prop = strrep(prop,'$2','ldata');
prop = strrep(prop,'$3','gdata');
prop = strrep(prop,'$4','ldata_time');
prop = strrep(prop,'$5','gdata_time');

%---------------------------------------------------------------------------
function umod = l_rparse_jac(umod)
%L_RPARSE_JAC Jacobian of propensities.

% "translate" into Matlab syntax:
prop = strrep(umod.private.rp.prop,'pow','power');
% (this is a bit of a hack since a known bug is that variables with
% names 'pow' or 'power' will ruin things...)
prop = strrep(umod.private.rp.prop,'!=','~=');

J = jacobian(str2sym(prop),str2sym(umod.private.rp.Species));
%[ii,jj,val] = find(J);
% the above gives unexpected behavior in some cases, the following
% seems to work better:
[ii,jj] = find(J ~= 0);
val = J(J ~= 0);
% for example, compare:
%syms a b;
%J = [a*b a == b];
%[ii,jj,val] = find(J)
%[ii,jj] = find(J ~= 0)
%val = J(J ~= 0)
ir = ii-1;
jc = cumsum([0 full(sparse(1,jj,1))]);
nnz = numel(val);
% again the following yields different results:
%ccode([a*b a == b])
%ccode(a == b)
% so seems better to do a "large" call to ccode
strval = ccode(val);
cval = strsplit(strval,'\n');
for i = 1:nnz
  % replace the destination of the result
  cval{i} = sprintf('val[%d] %s', ...
                    i-1,cval{i}(find(cval{i} == '=',1,'first'):end));
  % replace species and rates accordingly
  cval{i} = l_rewriteprop(cval{i},umod.private.rp.Species, ...
                          umod.private.rp.RateNames);
  cval{i} = l_finishprop(cval{i},umod.private.rp.Species, ...
                         umod.private.rp.RateNames, ...
                         umod.private.rp.ldrate,umod.private.rp.gdrate, ...
                         umod.private.rp.ldtrate,umod.private.rp.gdtrate);
end

% read file
[fid,msg] = fopen(umod.propensities,'rt');
if fid == -1, error(msg); end
F = char(fread(fid,64))';
if ~strncmp(F,'/* [Remove/modify this line not to overwrite this file] */',58)
  error('Will not overwrite existing file if not created by RPARSE.');
end
fclose(fid);

% construct a text-version Jacobian and write it
[fid,msg] = fopen(umod.propensities,'at');
if fid == -1, error(msg); end
F = '/* Sparse propensity Jacobian matrix. */\n';
F = [F ['/* Generated by RPARSE/L_RPARSE_JAC ' ...
        datestr(now,'yyyy-mm-dd HH:MM') ' */\n\n']];
F = [F sprintf('const int JAC_NNZ = %d; /* non-zeros in Jacobian */\n', ...
               nnz)];
F = [F sprintf('const int Mspecies = %d; /* number of species */\n\n', ...
               size(J,2))];
F = [F '/* static sparsity pattern */\n'];
F = [F 'static int JACir[] = {' sprintf('%d,',ir)];
F = [F(1:end-1) '};\nstatic int JACjc[] = {' sprintf('%d,',jc)];
F = [F(1:end-1) '};\n\n'];
F = [F '/* URDME sparse propensity Jacobian interface */\n'];
F = [F ...
     'void GET_jacobian(double *val,\n' ...
     '                  const URDMEstate_t *xstate,' ...
     'double time,double vol,\n' ...
     '                  const double *ldata,' ...
     'const double *gdata,\n' ...
     '                  const double *ldata_time,' ...
     'const double *gdata_time,int sd)\n{\n'];
for i = 1:nnz
  F = [F '  ' cval{i} '\n'];
end
F = [F '}\n'];
fprintf(fid,F);
fclose(fid);

%---------------------------------------------------------------------------
