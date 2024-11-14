%STATIONARY_BOUNDS Stationary analysis of Hes1.
%   This script analyses the stationary solution to Hes1@2 which,
%   since it is derived by (linear) stationary assumptions on Hes1@5,
%   gives results which are valid also in the latter case.
%
%   See also STATIONARY_STABILITY.

% S. Engblom 2024-07-09 (Revision)
% S. Engblom 2024-04-17

%% §1 estimates for N0 (= x0), the unique homogeneous stationary solution

% loop over range of parameters
k = 1; h = 4;
aa  = linspace(0.0,0.3,21)';
bb = [1e1 1e3 1e5];
x0 = zeros(numel(aa),numel(bb));
x1 = zeros(numel(aa),numel(bb),2);
opt = optimset('Display','off');
j = 0;
for b = bb
  j = j+1;
  i = 0;
  for a = aa'
    i = i+1;
    F = @(N)N.*(a+N.^k).*(1+b.*N.^h)-1;

    % upper bound for root
    % phi(N) > N*N^k*b*N^h-1 = 0
    uN(i,j) = b^(-1/(k+1+h));

    % lower bound for root
    % phi(N) < uN*(a+uN^k)*(1+b*N^h)-1 = 0
    lN(i,j) = ((1/uN(i,j)/(a+uN(i,j)^k)-1)/b)^(1/h);

    % find the root
    [x0(i,j),~,flag] = fzero(F,[lN(i,j) uN(i,j)],opt);
    if flag ~= 1, warning('Convergence issues.'); end
    
    % continue with the non-homogeneous roots
    f = @(x)1./(a+x.^k);
    g = @(x)1./(1+b*x.^h);
    F = @(N)[g(N(2))*f(N(1)); g(N(1))*f(N(2))]-N';
    [x1_,~,flag] = fsolve(F,[0.1 0.8],opt);
    if flag ~= 1, warning('Convergence issues.'); end
    x1(i,j,:) = reshape(x1_,1,1,2);
  end
end

% root
figure(1), clf,
plot(aa,x0,'b.-'); hold on,
title('Homogeneous steady-state');
xlabel('a');

% side test: expansion of the root x0 (see §7 at the end of this file)
eps_ = bb.^(-1/(k+h+1));
x0_ = arrayfun(@(a)polyval([(-a^2-1)/6 -a/3 1 0],eps_),aa, ...
               'UniformOutput',false);
x0_ = cat(1,x0_{:});
plot(aa,x0_,'b'); hold on,

% extra: upper and lower bounds of the root
% $$$ plot(aa,uN,'r');
% $$$ plot(aa,lN,'g');
% $$$ xlabel('a');
% $$$ for i = 1:numel(bb)
% $$$   text(aa(3),x0(3,i)-0.02,sprintf('b = %d',bb(i)));
% $$$ end

% homogeneous stationary solution unstable AND existence of 2-periodic
% stable solution provided this line is < 1:
f = @(x)1./(aa+x.^k);
g = @(x)1./(1+bb.*x.^h);
df = @(x)-k*x.^(k-1).*f(x).^2;
dg = @(x)-h*bb.*x.^(h-1).*g(x).^2;
% fg'-f'g < -1:
%bnd = f(x0).*dg(x0)-df(x0).*g(x0)+2;
% explicit dependence and somewhat nicer scaling:
bnd = h*bb.*x0.^h.\((1+bb.*x0.^h).*(1+k*x0.^(k+1).*(1+bb.*x0.^h)));
% (these relations are from Propositions 3.2 and 3.3)
plot(aa,bnd,'k--')
text(aa(end-10),bnd(end,1)-0.02,'< 1 implies non-homogeneous solution exists');
text(aa(2),bnd(end,1)+0.025,'< 1 implies homogeneous solution unstable');
axis([aa([1 end])' 0 1]);

% > 0 means non-homogeneous solution is stable
bnd = g(x1(:,:,1)).*g(x1(:,:,2)).*df(x1(:,:,1)).*df(x1(:,:,2))- ...
      g(x1(:,:,2)).*df(x1(:,:,1))- ...
      g(x1(:,:,1)).*df(x1(:,:,2))- ...
      dg(x1(:,:,1)).*dg(x1(:,:,2)).*f(x1(:,:,1)).*f(x1(:,:,2))+1;

%% §2 "pedagogical" figure for some suitable parameters
figure(2), clf, htl = tiledlayout(1,2); nexttile,
a = 1;
b = 10;
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b*x.^h);
phi = @(x)f(x).*g(x);

x = linspace(0,1);
plot(x,phi(x),'r','linewidth',2); hold on,
plot(x,x,'b','HandleVisibility','off');
opt = optimset('Display','off','tolfun',1e-12,'tolx',1e-12);
F = @(N)phi(N)-N;
x0 = fzero(F,[0 1],opt);
plot(x0,x0,'ro','HandleVisibility','off');
legend('$\varphi$','Interpreter','latex');

% make this plot publishable
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',[0 1],'ytick',[1],'fontsize',9);

%% §3 existence of 2-periodic solution (2-cell problem)
% 0 = D2f(N1)-N1 } N1/f(N1) = g(N2)
% 0 = g(N1)-D1   }
% 0 = D1f(N2)-N2 } N2/f(N2) = g(N1)
% 0 = g(N2)-D2   }

% parameters a bit shifted for visual purposes:
a = 0.1;
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b*x.^h);

F = @(N)N.*(a+N.^k).*(1+b.*N.^h)-1;
x0 = fzero(F,[0 1],opt);

F = @(N)[g(N(2))*f(N(1)); g(N(1))*f(N(2))]-N';
x1 = fsolve(F,[0.1 0.8],opt);

% continue "pedagogical" figure, plotting to the right
nexttile,
plot(x,x,'b','HandleVisibility','off'); hold on,

% h(x) := x/f(x) = g(x) or x = invh(g(x)) =: gam1(x), search for root
% x in [0,1] of h(x) = g(x):
invh_ = @(y)fzero(@(x)x/f(x)-y,[0 1],opt); % scalar
gam1 = @(x)arrayfun(invh_,g(x)); % vectorized
gam2 = @(x)gam1(gam1(x));

plot(x,gam2(x),'b','LineWidth',2);
plot(x1,gam2(x1),'bo','HandleVisibility','off');
plot(x0,gam2(x0),'ro','HandleVisibility','off');
plot(x,gam1(x),'r--','LineWidth',2);
legend('$\gamma_2$','$\gamma$','Location','SE', ...
       'Interpreter','latex');

% make this plot publishable
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',[0 1],'ytick',[],'fontsize',9);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[200 200 400 120]);
htl.Padding = 'compact';
htl.TileSpacing = 'compact';
%print -depsc ~/Desktop/fixpnt.eps

return;

%% §4 solving the limiting inequality directly
k = 1; h = 4;
X0 = @(a,b)fzero(@(N)N.*(a+N.^k).*(1+b.*N.^h)-1,[0 1],opt);
B0 = @(a)fsolve(@(b)(1+b*X0(a,b)^h)*(1+k*X0(a,b)^(k+1)*(1+b*X0(a,b)^h))- ...
                b*h*X0(a,b)^h,4.0,opt);
% b > B0(a) implies 2-periodic solution exists
a = linspace(0,2,21);
b0 = arrayfun(B0,a);
figure(3), clf,
plot(a,b0);
xlabel('a'); ylabel('b_0');
title('if b > b0(a), then non-homogeneous steady-state exists');

return;

%% §5 phase portrait in the (N1,N2)-plane
N1 = linspace(0,1,41);
N2 = linspace(0,1,39);
[N1_,N2_] = meshgrid(N1,N2);
% stationary in (D1,D2)
D1 = g(N1_);
D2 = g(N2_);
F1 = D2.*f(N1_)-N1_;
F2 = D1.*f(N2_)-N2_;

% normalize
d = sqrt(F1.^2+F2.^2);
F1 = F1./d;
F2 = F2./d;

figure(4), clf,
h = quiver(N1_,N2_,F1,F2,0.6);
hold on,
plot(x1,x1(end:-1:1),'bo','HandleVisibility','off');
plot(x0,x0,'ro','HandleVisibility','off');
axis([0 1 0 1]);

return;

%% §6 attempt to find more than one cyclic solution by varying (a,b,k,h)
figure(5); clf,
h1 = []; h2 = [];
k = 3; h = 5; a = 1; b = 1e2;

while 1
  f = @(x)1./(a+x.^k);
  g = @(x)1./(1+b*x.^h);
  invh_ = @(y)fzero(@(x)x/f(x)-y,[0 1],opt);
  gam1 = @(x)arrayfun(invh_,g(x));
  gam2 = @(x)gam1(gam1(x));

  delete(h1); delete(h2);
  plot(x,x,'k'); hold on,
  h1 = plot(x,gam2(x),'b');
  h2 = plot(x,gam1(x),'r--');
  key = input('a/A/b/B: ','s');
  switch key,
    case 'a', a = a/2
    case 'A', a = a*2
    case 'b', b = b/1.25
    case 'B', b = b*1.25
    otherwise, break;
  end
end
% (seems not doable!)

return;

%% §7 analytic expansion of x0 in terms of a and b

% 1 = x*(a+x^k)*(1+b*x^h)
%
% b*x^(k+h+1) = 1-a*x*(1+b*x^h)-x*(a+x^k)
syms a b eps_ x x_ c1 c2 c3 c4 c5

% solve for steady-state solutions as an expansion
k = 1; h = 4;
assume(eps_ > 0); % eps_ = b^(-1/(k+h+1));

g1 = @(x)b*x^(k+h+1);
g2 = @(x)1-a*x*(1+b*x^h)-x*(a+x^k);

n = 5;
x_ = 0;
for j = 1:n
  for i = 1:j
    x_ = series(simplify(eps_*g2(x_)^(1/(k+h+1))),eps_,'order',j);
  end
  c = simplify(coeffs(x_,eps_));
end

% conclusion: x0 ~ polyval([(-a^2-1)/6 -a/3 1 0],eps_) with eps_ =
% b^(-1/(k+h+1)) is considered small
