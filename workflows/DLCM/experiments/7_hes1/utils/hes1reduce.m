function [a,b,v,tau,M0_,D0_] = hes1reduce(alpha,mu,H)
%HES1REDUCE Reduced Hes1-model given full parameters.
%   [a,b,v,tau,M0_,D0_] = HES1REDUCE(alpha,mu,H) essentially does the
%   same job as HES1RED_PARAMS for alternative == 5. The intent with
%   this function is to work as a convenient auxiliary function in
%   various other functions.
%
%   No error-checking is performed.
%
%   See also HES1UNREDUCE, HES1_PARAMS, HES1RED_PARAMS.

% S. Engblom 2024-10-28

c1 = alpha(2)*alpha(3)/mu(2);
c2 = (alpha(4)/(mu(4)*H(1)))^H(3);
c3 = alpha(1)*alpha(5)/mu(5);
c4 = (alpha(4)/(mu(4)*H(2)))^H(4);
D0_ = c3/mu(1);
M0_ = (c1*D0_/(mu(3)*c2))^(1/(H(3)+1));

a = 1/(c2*M0_^H(3));
b = c4*M0_^H(4);
v = mu(1)/mu(3);
tau = 1/mu(3);
