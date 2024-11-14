function [x0,x12] = hes1red_nonhomogeneous(a,b,k,h)
%HES1RED_NONHOMOGENEOUS Non-homogeneous solution.
%   [X0,X12] = HES1RED_NONHOMOGENEOUS(a,b,k,h) returns both the
%   homogeneous solution X0 and the non-homogeneous solution X12 = [X1
%   X2] for reduced parameters (a,b,k,h).
%
%   This function solves for the stable non-homogeneous solution
%   only.
%
%   See also HES1RED_HOMOGENEOUS.

% S. Engblom 2024-09-30

% solve for homogeneous solution first
phi = @(x)x.*(a+x.^k).*(1+b.*x.^h)-1;
opts = optimset('Display','off','tolx',1e-12,'tolfun',1e-12);
[x0,~,flag] = fzero(phi,[0 1],opts);
if flag ~= 1
  disp(sprintf('Convergence issuses, flag = %d.',flag));
end

% next define the more difficult non-homogeneous problem
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b.*x.^h);
invH = @(y)fzero(@(x)x/f(x)-0.5*g(x)-y,[0 1],opts);
Gam1 = @(x)invH(0.5*g(x));
invh = @(y)fzero(@(x)x/f(x)-y,[0 1],opts);
gam1 = @(x)invh(g(x));
%Gam12 = @(x)Gam1(gam1(x));
Gam21 = @(x)gam1(Gam1(x));

% bracket the roots & solve
xx = [linspace(0,(1-1e-4)*x0,50)' linspace((1+1e-4)*x0,1,50)'];
[ii,jj] = find(diff(sign(arrayfun(Gam21,xx)-xx)));
ij = sub2ind(size(xx),ii,jj);
if numel(ij) == 0
  x12 = nan(1,2);
  return;
end
[x_,~,flag] = fzero(@(x)Gam21(x)-x, ...
                    xx([ij(1) ij(1)+1]),opts);
if flag ~= 1
  disp(sprintf('Convergence issuses, flag = %d.',flag));
  x12 = nan(1,2);
  return;
end

x12 = [x_ Gam1(x_)];
