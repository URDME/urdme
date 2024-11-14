function J = hes1_Jacobian(X,alpha,mu,H,funJ,funC,W)
%HES1_JACOBIAN Jacobian for the full Hes1-model.
%   HES1_JACOBIAN works exactly like HES1_SYSTEM but returns the
%   Jacobian matrix instead. Also, note that the input "fun" should
%   rather be "funJ", see HES1_BUILDODE.
%
%   See also HES1_BUILDODE, HES1_SYSTEM.

% S. Engblom 2024-09-30

% by the linearity of the couplings, they can be represented by a
% matrix, fully known once the parameters have been proposed:
C = kron(W,funC(alpha,mu,H));

% assemble
X = reshape(X,[],size(W,2));
[m,n] = size(X);
J = zeros(0,m*n);
for i = 1:n
  J = [J; zeros(m,(i-1)*m) funJ(X(:,i),alpha,mu,H) zeros(m,(n-i)*m)];
end
J = J+C;
