function [expr,data] = seqexpand(expr,data,seq)
%SEQEXPAND Expand sequences of expressions.
%   EXPR = SEQEXPAND(EXPR,SEQ) expands the expressions EXPR according
%   to the sequence SEQ. EXPR is a cell-vector of strings (or a single
%   string) and SEQ is a cell-vector of alternatively sequence
%   variable names (single characters) and numerical sequences (double
%   vectors).
%
%   The resulting cell-vector EXPR is the result of replacing all
%   occurences of '$<sequence variable>' with all sequence values,
%   according to the ordering in SEQ. For example,
%   SEQEXPAND('($i)*($j)',{'i' 1:2 'j' [-1 1]}) returns {'(1)*(-1)'
%   '(2)*(-1)' '(1)*(1)' '(2)*(1)'}.
%
%   Dependent variables may be expressed by stacking them on top of
%   each other, for example, SEQEXPAND('K$i$j',{['ij']' [1:4;
%   4:-1:1]}) returns {'K14' 'K23' 'K32' 'K41'}.
%
%   [EXPR,DATA] = SEQEXPAND(EXPR,DATA,SEQ) does the same thing but
%   simultaneously distributes a corresponding cell-vector DATA
%   containing numerical arrays. Scalars are distributed by making
%   copies and arrays by dividing them evenly over the resulting
%   expressions. For example, SEQEXPAND('K$i$j',1:4,{'i' 1:2}) will
%   distribute the data into the chunks [1 2] and [3 4], respectively.
%
%   Note: this is a utility function and limited error-checking is
%   performed.
%
%   Examples:
%     % mutiple sequences, including a dependent variable 'e'
%     expr = seqexpand('x$i + y$e * z$i$j', ...
%                      {['ie']' [1:3; 3:-1:1] 'j' 1:5});
%
%     % distribution of data
%     [expr,data] = seqexpand({'x$i$j' 'y$i$j' 'z$i'},{1:6 1:2 1:6}, ...
%                             {'i' 1:2 'j' 1:3});

% S. Engblom 2017-04-14

% syntax
if nargin < 3
  seq = data;
  data = [];
end

% expression
if ~iscellstr(expr)
  if ~ischar(expr)
    error('Expression should be strings.');
  end
  expr = {expr};
end

% data, if any
if ~isempty(data)
  if ~iscell(data)
    data = {data};
  end
  if any(size(expr) ~= size(data))
    error('Expressions and data must match in size,');
  end
end

% sequence representation
vars = seq(1:2:end);
if ~iscellstr(vars) || any(cellfun('size',vars,2) ~= 1)
  error('All sequence variable names must be single characters.');
end
seqs = seq(2:2:end);
lens = cellfun('size',seqs,2);
  
% loop over expr, expanding each expr{i} in turn
for i = 1:numel(expr)
  % expand "in place" using an expanding cell vector
  expr{i} = {expr{i}};
  if ~isempty(data)
    data{i} = {data{i}(:)};
  end

  % loop over all sequence variables in seqs
  for j = 1:numel(vars)
    % if sequence variable was found, expand it
    if ~isempty(intersect(expr{i}{1}(find(expr{i}{1} == '$')+1),vars{j}))
      % expand expr{i} to the correct size
      expr{i} = repmat(expr{i},1,lens(j));

      % divide data among the resulting expressions
      if ~isempty(data)
        if isscalar(data{i}{1})
          data{i} = reshape(repmat(data{i},1,lens(j)),[],1);
        else
          siz = numel(data{i}{1})/lens(j);
          if siz ~= fix(siz)
            error('Data can not be evenly distributed.');
          end
          siz = repmat(siz,1,lens(j));
          % distribute all data using cell vectors, then concatenate to one cell
          data{i} = cellfun(@(x)(mat2cell(x,siz)),data{i},'uniformoutput',false);
          data{i} = cat(1,data{i}{:});
        end
      end

      % loop over all dependent sequence variables
      for k = 1:size(vars{j},1)
        % distribute sequence variables using a cell array, then convert each
        % to strings and replace all at once
        replace = repmat(num2cell(seqs{j}(k,:)),size(expr{i},1),1);
        replace = cellfun(@(x){num2str(x)},replace);
        expr{i} = strrep(expr{i},['$' vars{j}(k)],replace);
      end
      expr{i} = reshape(expr{i},[],1);
    end
  end
end

% "flatten out" to one row list of string expressions
expr = cat(1,expr{:})';
if ~isempty(data)
  data = cat(1,data{:})';
end
