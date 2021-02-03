function varargout = map(L,G,varargin)
%MAP Map indices.
%   [L1,L2,...] = MAP(L,G,G1,G2,...) maps the global indices G1, G2,
%   ... into local indices L1, L2, ..., following the map G --> L.
%
%   All index vectors are assumed to be column vectors.
%
%   Example:
%     g = [5 5 5 3 3 1]';
%     l = map([1 2 3 4 5]',[5 4 3 2 1]',g)

% S. Engblom 2017-12-20 (minor revision)
% S. Engblom 2017-10-06

% hashing used to find Gi in G
Gi = cat(1,varargin{:});
[~,ix] = fsetop('ismember',Gi',G');

% then perform the map G --> L
str = cumsum([1 cellfun('prodofsize',varargin)]);
varargout = cell(1,numel(varargin));
for i = 1:numel(varargout)
    varargout{i} = L(ix(str(i):str(i+1)-1));
end
