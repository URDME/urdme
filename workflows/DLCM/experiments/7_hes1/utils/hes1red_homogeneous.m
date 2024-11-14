function x0 = hes1red_homogeneous(a,b,k,h)
%HES1RED_HOMOGENEOUS Homogeneous solution.
%   X0 = HES1RED_HOMOGENEOUS(a,b,k,h) solves for the unique
%   homogeneous solution of the reduced Hes1-model. Inputs are the
%   reduced parameters (a,b,k,h) and the function returns X0. The
%   other variable is obtained as Y0 = 1/(1+b*X0^h).
%
%   See also HES1RED_NONHOMOGENEOUS.

% S. Engblom 2024-09-30

% solve
phi = @(x)x.*(a+x.^k).*(1+b.*x.^h)-1;
opts = optimset('Display','off','tolx',1e-12,'tolfun',1e-12);
[x0,~,flag] = fzero(phi,[0 1],opts);
if flag ~= 1
  disp(sprintf('Convergence issuses, flag = %d.',flag));
end
