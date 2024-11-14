function X = hes1unreduce(x,y,M0_,D0_,alpha,mu,H,W)
%HES1UNREDUCE Map to full Hes1-model from reduced case.
%   X = HES1UNREDUCE(x,y,M0_,D0_,alpha,mu,H,W) maps to states X = [D N
%   M P n] from states (x,y) using the parameters
%   (M0_,D0_,alpha,mu,H). Also, W is the coupling matrix for the
%   considered model. The intent with this function is to work as a
%   convenient auxiliary function in various other functions.
%
%   No error-checking is performed.
%
%   See also HES1REDUCE, HES1_PARAMS, HES1RED_PARAMS.

% S. Engblom 2024-10-28

M = x(:)'*M0_;
D = y(:)'*D0_;
N = alpha(2)/mu(2)*D*W';
P = alpha(4)/mu(4)*M;
n = alpha(5)/mu(5)./(1+(P/H(2)).^H(4));
X = [D; N; M; P; n];
