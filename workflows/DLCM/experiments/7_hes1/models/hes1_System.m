function F = hes1_System(X,alpha,mu,H,fun,funC,W)
%HES1_SYSTEM RHS function for the Hes1 system on a grid.
%   F = HES1_SYSTEM(X,alpha,mu,H,fun,funC,W) returns the RHS F for
%   states X of the Hes1-model. Parameters are alpha = alpha_[D N M M
%   n], similarly for mu and for H = [KM Kn k h]. The RHS is
%   preassembled by fun and funC as obtained from HES_BUILDODE and W
%   is the coupling matrix over the grid, usually satisfying sum(W,2)
%   = 1.
%
%   See also HES1_BUILDODE, HES1_JACOBIAN.

% S. Engblom 2024-09-30

% by the linearity of the couplings, they can be represented by a
% matrix, fully known once the parameters have been proposed:
C = kron(W,funC(alpha,mu,H));

% RHS evaluation by independent local evaluations...
X = reshape(X,[],size(W,2));
F = zeros(size(X));
for i = 1:size(X,2)
  F(:,i) = fun(X(:,i),alpha,mu,H);
end

% ...followed by the linear coupling
F = F(:)+C*X(:);
