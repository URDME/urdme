function G = G_inline(K,I,N)
%G_INLINE Determine dependency graph for inline propensities.
%   G = G_INLINE(K,I,N) returns the dependency graph G for the inline
%   propensities defined by (K,I,N). No error-checking is performed.
%
%   See NSM.

% S. Engblom 2017-03-21

% immediate
i = (K([1 1 2],:) ~= 0).*I;
j = repmat(1:size(N,2),3,1);
keep = find(i ~= 0);
i = i(keep);
j = j(keep);
H = sparse(i,j,1,size(N,1),size(N,2));

% [diffusion reaction]-parts of G
G = [H' H'*abs(N)];
G = double(G ~= 0);
