%STATIONARY_3CELL Stationarity of Hes1-model for 3 cells.
%   This script investigates the stability of the reduced Hes1-model
%   over three cells.

% S. Engblom 2024-09-11 (stationary_3cell)
% S. Engblom 2024-08-22 (stationary_graph)
% S. Engblom 2024-03-12
% S. Engblom 2024-03-02 (hes1_stability on a grid)
% S. Engblom 2024-01-20 (hes1_stability)

%% (1) analysis of the 2- or 3-cell problems
% $$$ n = 2;
% $$$ syms W [n n]
% $$$ if n == 2
% $$$   W = [0 1; 1/2 1/2];
% $$$ else
% $$$   W = [0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0];
% $$$ end
% $$$ 
% $$$ % the model
% $$$ syms f(x) g(x) df(x) dg(x)
% $$$ syms x [n 1]
% $$$ F = (W*g(x)).*f(x)-x;
% $$$ 
% $$$ % Jacobian
% $$$ J = diag(f(x))*W*diag(dg(x))+diag(W*g(x).*df(x))-eye(n);
% $$$ % check of this:
% $$$ % $$$ J0 = jacobian(F,x);
% $$$ % $$$ for i = 1:n
% $$$ % $$$   J0 = subs(J0,[diff(f(x(i)),x(i))],[df(x(i))]);
% $$$ % $$$   J0 = subs(J0,[diff(g(x(i)),x(i))],[dg(x(i))]);
% $$$ % $$$ end
% $$$ % $$$ J-J0 % = 0
% $$$ 
% $$$ % stability of homogeneous and non-homogeneous solutions
% $$$ syms x0
% $$$ lam0 = eig(subs(J,x,repmat(x0,n,1)))
% $$$ % unfortunately, these are not very telling:
% $$$ % $$$ if n == 2
% $$$ % $$$   lam12 = eig(J);
% $$$ % $$$ else
% $$$ % $$$   lam123 = eig(subs(J,x,[x(1) x0 x(2)].'))
% $$$ % $$$   lam112 = eig(subs(J,x,x([1 1 2])))
% $$$ % $$$   lam122 = eig(subs(J,x,x([1 2 2])))
% $$$ % $$$ 
% $$$ % $$$   syms a b h k
% $$$ % $$$   lam112 = subs(lam112,dg(x),-b*h*x.^(h-1).*g(x).^2);
% $$$ % $$$   lam112 = subs(lam112,df(x),-k*x.^(k-1).*f(x).^2);
% $$$ % $$$   syms invW
% $$$ % $$$   invW = inv(W);
% $$$ % $$$   lam112 = subs(lam112,g(x),invW*(x./f(x)));
% $$$ % $$$   lam112 = subs(lam112,x,x([1 1 2]))
% $$$ % $$$ 
% $$$ % $$$ % det = prod(n eigenvalues)
% $$$ % $$$ dJ = det(J);
% $$$ % $$$ simplify(subs(dJ,x,repmat(x0,n,1)))
% $$$ % $$$ simplify(subs(dJ,x,x([1 2 2])))
% $$$ % $$$ end

%% (2) 3-cell problem
if ~exist('save_data','var')
  save_data = false;
end

par = hes1_params;
alpha = [par.alphaD par.alphaN par.alphaM par.alphaP par.alphan];
mu = [par.muD par.muN par.muM par.muP par.mun];
H = [par.KM par.Kn par.k par.h];
k = par.k;
h = par.h;
W = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
opts = optimset('Display','off','tolx',1e-8,'tolfun',1e-10);

% determine critical bifurcation scaling
crit = @(s)l_criterion( ...
    [par.alphaD par.alphaN/s par.alphaM par.alphaP par.alphan], ...
    [par.muD par.muN par.muM par.muP par.mun],[par.KM par.Kn par.k par.h]);
opts = optimset('Display','off');
[scrit,~,flag] = fzero(crit,50,opts);
if flag ~= 1, warning('Convergence issues.'); end

crit2 = @(s)l_criterion2( ...
    [par.alphaD par.alphaN/s par.alphaM par.alphaP par.alphan], ...
    [par.muD par.muN par.muM par.muP par.mun],[par.KM par.Kn par.k par.h]);
[scrit2,~,flag] = fzero(crit2,[scrit 50],opts);
if flag ~= 1, warning('Convergence issues.'); end

% sweep
svec = linspace(1,100);
svec = [svec(svec < scrit) scrit svec(scrit < svec & svec < scrit2) ...
        scrit2 svec(scrit2 < svec)];
x0 = nan(1,numel(svec));
x1 = nan(2,numel(svec));
x2 = nan(2,numel(svec));
% saved state "P" in full model:
X0 = nan(1,numel(svec));
X1 = nan(2,numel(svec));
X2 = nan(2,numel(svec));
lam = nan(3,numel(svec));
res = nan(3,numel(svec));
% solutions ordered according to:
% x0 x1a x1b
% x0 x2a x2b
% x0 x2a x2b
for i = 1:numel(svec)
  s = svec(i);

  % reduced parameters (a,b)
  [a,b,~,~,M0_,D0_] = hes1reduce(alpha./[1 s 1 1 1],mu,H);

  % useful defines
  f = @(x)1./(a+x.^k);
  g = @(x)1./(1+b*x.^h);
  df = @(x)-k*x.^(k-1)./(a+x.^k).^2;
  dg = @(x)-b*h*x.^(h-1)./(1+b*x.^h).^2;
  F = @(x)(W*g(x)).*f(x)-x;
  J = @(x)diag(f(x))*W*diag(dg(x))+diag(W*g(x).*df(x))-eye(3);

  % (1) homogeneous solution
  x0(i) = hes1red_homogeneous(a,b,k,h);

  % next define the more difficult non-homogeneous problem
  invH = @(y)fzero(@(x)x/f(x)-0.5*g(x)-y,[0 1],opts);
  Gam1 = @(x)invH(0.5*g(x));
  invh = @(y)fzero(@(x)x/f(x)-y,[0 1],opts);
  gam1 = @(x)invh(g(x));
  % 0 = g(x2)f(x1)-x1
  % 0 = (0.5*g(x1)+0.5*g(x2))f(x2)-x2
  %
  % x1/f(x1) =: h(x1) = g(x2)
  % x2/f(x2)-0.5*g(x2) =: H(x2) = 0.5*g(x1)
  %
  % x1 = invh(g(x2)) = gam1(x2)
  % x2 = invH(d*g(x1)) = Gam1(x1)
  %
  % x2 = Gam1(gam1(x2)) =: Gam12(x2)
  % x1 = gam1(Gam1(x1)) =: Gam21(x1)
  Gam12 = @(x)Gam1(gam1(x));
  Gam21 = @(x)gam1(Gam1(x));

  % (2) non-homogeneous solution
  xx = [linspace(0,(1-1e-4)*x0(i),50)' linspace((1+1e-4)*x0(i),1,50)'];
  [ii,jj] = find(diff(sign(arrayfun(Gam21,xx)-xx)));
  ij = sub2ind(size(xx),ii,jj);
  for j = 1:numel(ij)
    [x_,~,flag] = fzero(@(x)Gam21(x)-x, ...
                        xx([ij(j) ij(j)+1]),opts);
    if flag ~= 1
      disp(sprintf('Convergence issuses, flag = %d.',flag));
    end
    x1(j,i) = x_;
  end
  x2(1:numel(ij),i) = arrayfun(Gam1,x1(1:numel(ij),i));

  % (2) spectra
  xroot = [x0(i) x0(i) x0(i)]';
  X_ = hes1unreduce(xroot,g(xroot),M0_,D0_,alpha./[1 s 1 1 1],mu,H,W);
  X0(i) = X_(4,1);
  J_ = J(xroot);
  lam(1,i) = max(real(eig(J_)));
  res(1,i) = norm(F(xroot));
  xroot = [x1(1,i) x2(1,i) x2(1,i)]';
  X_ = hes1unreduce(xroot,g(xroot),M0_,D0_,alpha./[1 s 1 1 1],mu,H,W);
  X1(:,i) = X_(4,1:2)';
  J_ = J(xroot);
  lam(2,i) = max(real(eig(J_)));
  res(2,i) = norm(F(xroot));
  xroot = [x1(2,i) x2(2,i) x2(2,i)]';
  X_ = hes1unreduce(xroot,g(xroot),M0_,D0_,alpha./[1 s 1 1 1],mu,H,W);
  X2(:,i) = X_(4,1:2)';
  J_ = J(xroot);
  lam(3,i) = max(real(eig(J_)));
  res(3,i) = norm(F(xroot));
end

% $$$ figure(5), clf, hold on,
% $$$ ix = find(lam(2,:) > 0);
% $$$ plot(svec(ix),[x1(1,ix); x2(1,ix)],'b--', ...
% $$$      'LineWidth',2,'HandleVisibility','off');
% $$$ ix = find(lam(2,:) <= 0);
% $$$ plot(svec(ix),x1(1,ix),'b', ...
% $$$      'LineWidth',2,'HandleVisibility','off');
% $$$ plot(svec(ix),x2(1,ix),'b', ...
% $$$      'LineWidth',2);
% $$$ 
% $$$ ix = find(lam(3,:) > 0);
% $$$ plot(svec(ix),[x1(2,ix); x2(2,ix)],'b--', ...
% $$$      'LineWidth',2,'HandleVisibility','off');
% $$$ ix = find(lam(3,:) <= 0);
% $$$ plot(svec(ix),[x1(2,ix); x2(2,ix)],'b', ...
% $$$      'LineWidth',2,'HandleVisibility','off');
% $$$ 
% $$$ ix = find(lam(1,:) > 0);
% $$$ plot(svec(ix),x0(ix),'r--','LineWidth',2, ...
% $$$      'HandleVisibility','off')
% $$$ ix = find(lam(1,:) <= 0);
% $$$ plot(svec(ix),x0(ix),'r','LineWidth',2);
% $$$ xline(scrit,'k','HandleVisibility','off');
% $$$ xline(scrit2,'k-.','HandleVisibility','off');

% as above, but plot in terms of Hes1@5 states
figure(1), clf,
htl = tiledlayout(1,2); nexttile,
% scale to [0,1]:
Pmax = max(X0(:),[],'omitnan');
Pmax = max(Pmax,max(X1(:),[],'omitnan'));
Pmax = max(Pmax,max(X2(:),[],'omitnan'));
X0 = X0/Pmax;
X1 = X1/Pmax;
X2 = X2/Pmax;
ix = find(lam(2,:) > 0);
plot(svec(ix),[X1(1,ix); X1(2,ix)],'b--', ...
     'LineWidth',2,'HandleVisibility','off');
hold on,
set(gca,'YScale','log'); % for consistenty with other plots
ix = find(lam(2,:) <= 0);
plot(svec(ix),X1(1,ix),'b', ...
     'LineWidth',2,'HandleVisibility','off');
plot(svec(ix),X1(2,ix),'b', ...
     'LineWidth',2);

ix = find(lam(3,:) > 0);
plot(svec(ix),[X2(1,ix); X2(2,ix)],'b--', ...
     'LineWidth',2,'HandleVisibility','off');
ix = find(lam(3,:) <= 0);
plot(svec(ix),[X2(1,ix); X2(2,ix)],'b', ...
     'LineWidth',2,'HandleVisibility','off');

ix = find(lam(1,:) > 0);
plot(svec(ix),X0(ix),'r--','LineWidth',2, ...
     'HandleVisibility','off')
ix = find(lam(1,:) <= 0);
plot(svec(ix),X0(ix),'r','LineWidth',2);
xline(scrit,'k','HandleVisibility','off');
xline(scrit2,'k-.','HandleVisibility','off');

% make this plot publishable
legend('Non-homogeneous','Homogeneous', ...
       'Interpreter','latex','location','SE');
set(gcf,'defaultTextInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');
xlabel('$s$');
ylabel('$P$');
axis([svec(1) svec(end) 1e-3 1]);
set(gca,'xtick',0:50:svec(end),'ytick',0:1,'yticklabel',1,'fontsize',9);
box on

if save_data
  disp('Saving...')
  save('../data/stationary_3cell.mat', ...
       'svec','scrit','scrit2','x0','x1','x2','lam', ...
       'X0','X1','X2');
  disp('   ...saved.')
end

%% (3) pedagogical figure illustrating the condition for existence
nexttile,

% reload parameters (avoids saving them in stationary_3cell.mat):
par = hes1_params;
alpha = [par.alphaD par.alphaN par.alphaM par.alphaP par.alphan];
mu = [par.muD par.muN par.muM par.muP par.mun];
H = [par.KM par.Kn par.k par.h];
k = par.k;
h = par.h;
W = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
opts = optimset('Display','off','tolx',1e-8,'tolfun',1e-10);

ss = [47 100];
for j = 1:numel(ss)
  [~,i] = min(abs(ss(j)-svec));
  s = svec(i);

  % shift the appearance to center on the homogeneous solution:
  if j == 1
    x0_ = x0(i);
    xx = linspace(0.1,0.6,101);
    xx = [xx(xx < x0_) x0_ xx(xx > x0_)];
    plot(xx,xx,'b','HandleVisibility','off'); hold on,
    xl = xx(xx <= x0_);
    xr = xx(xx >= x0_);
    dx = 0;
  else
    dx = x0_-x0(i);
  end

  [a,b] = hes1reduce(alpha./[1 s 1 1 1],mu,H);
  f = @(x)1./(a+x.^k);
  g = @(x)1./(1+b*x.^h);
  invH = @(y)fzero(@(x)x/f(x)-0.5*g(x)-y,[0 1],opts);
  Gam1 = @(x)invH(0.5*g(x));
  invh = @(y)fzero(@(x)x/f(x)-y,[0 1],opts);
  gam1 = @(x)invh(g(x));
  Gam12 = @(x)Gam1(gam1(x));
  Gam21 = @(x)gam1(Gam1(x));

  if j == 1
    % increase the visibility by scaling to make the zeros clearer: 
    xx = [xl xr];
    yy = [arrayfun(Gam21,xl) arrayfun(Gam12,xr)];
    % rotate the coordinate system...
    theta = pi/4;
    R = [cos(theta) -sin(theta); ...
         sin(theta) cos(theta)];
    zz = R*[xx; yy];
    % ...scale in the appropriate direction...
    zz(1,:) = 20*zz(1,:);
    % ...and rotate back:
    zz = R'*zz;
    plot(zz(1,:),zz(2,:),'b','LineWidth',2);
    % unscaled:
    %plot(xl,arrayfun(Gam21,xl),'b','LineWidth',2);
    %plot(xr,arrayfun(Gam12,xr),'b','LineWidth',2,'HandleVisibility','off');
    plot(x0(i),x0(i),'ro','HandleVisibility','off');
    plot([x1(1,i) x2(1,i)],[x1(1,i) x2(1,i)],'bo','HandleVisibility','off');
    plot([x1(2,i) x2(2,i)],[x1(2,i) x2(2,i)],'k+','HandleVisibility','off');
  elseif j == 2
    plot(xl,arrayfun(Gam21,xl-dx)+dx,'b','LineWidth',2,'HandleVisibility','off');
    plot(xr,arrayfun(Gam12,xr-dx)+dx,'b','LineWidth',2,'HandleVisibility','off');
  end
end

% make this plot publishable
axis([0.2 0.52 0.2 0.52]);
legend('$\Gamma_2$','Location','NW', ...
       'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',[],'ytick',[],'fontsize',9);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[200 200 400 120]);
htl.Padding = 'compact';
htl.TileSpacing = 'compact';
%print -depsc ~/Desktop/3cell_bifurcation.eps

%% 3-cell phase portrait in a selected plane
i = 1;
s = svec(i);

% these 3 points select the plane:
p0 = [x0(i) x0(i) x0(i)];
p1 = [x1(1,i) x2(1,i) x2(1,i)];
p2 = [x1(2,i) x2(2,i) x2(2,i)];
v1 = p1-p0;
v2 = p2-p0;
% so the plane is (x,y,z) = p0+a*v1+b*v2

% reduced parameters (a,b)
[a,b,v,tau,M0_,D0_] = hes1reduce(alpha./[1 s 1 1 1],mu,H);
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b*x.^h);
  
X1 = linspace(0,1,15);
X2 = linspace(0,1,17)';
[X1_,X2_] = meshgrid(X1,X2);

% find (a,b) from the given (x,y)-coordinates in the plane:
ab = [v1(1:2)' v2(1:2)']\[X1_(:)-p0(1) X2_(:)-p0(2)]';
% then use this to deduce the z-coordinate:
X3_ = reshape(p0(3)+ab(1,:)*v1(3)+ab(2,:)*v2(3),size(X1_));

% compute RHS
Y1 = g(X1);
Y2 = g(X2);
Y3 = g(X3_);
F1 = (0.5*Y2+0.5*Y3).*f(X1)-X1;
F2 = (0.5*Y1+0.5*Y3).*f(X2)-X2;
F3 = (0.5*Y1+0.5*Y2).*f(X3_)-X3_;

figure(2), clf,
quiver(X1_,X2_,F1,F2,3);
hold on,

% a finer resolution for the contours
X1 = linspace(0,1,63);
X2 = linspace(0,1,64)';
[X1_,X2_] = meshgrid(X1,X2);
ab = [v1(1:2)' v2(1:2)']\[X1_(:)-p0(1) X2_(:)-p0(2)]';
X3_ = reshape(p0(3)+ab(1,:)*v1(3)+ab(2,:)*v2(3),size(X1_));
Y1 = g(X1);
Y2 = g(X2);
Y3 = g(X3_);
F1 = (0.5*Y2+0.5*Y3).*f(X1)-X1;
F2 = (0.5*Y1+0.5*Y3).*f(X2)-X2;
F3 = (0.5*Y1+0.5*Y2).*f(X3_)-X3_;

% magnitude of RHS as contours
D = sqrt(F1.^2+F2.^2+F3.^2);
contour(X1_,X2_,D,2.^-[0 0.5 1.5 2.25],'b')
axis([0 1 0 1]);

% add sample trajectories
[fun,funJ,funC] = hes1_buildODE;
Tend = 24*60*1.5;
tspan = linspace(0,Tend);

% build initial data
Y0_ = [0.95 0.975 1; 0.9 0.95 0.975]; % ~upper right corner

for j = 1:2
  Y0 = Y0_(j,:);
  Y0 = [Y0; g(Y0)]; % (x,y) in 3-cell system
  % go to full system:
  Y0 = hes1unreduce(Y0(1,:),Y0(2,:),M0_,D0_,alpha./[1 s 1 1 1],mu,H,W);

  % solve
  [~,Y] = ode15s(@(t,y)hes1_System(y,alpha,mu,H, ...
                                   fun,funC,W),tspan,Y0, ...
                 odeset('AbsTol',1e-8,'RelTol',1e-10));
  Y = Y';

  % truncate back to the plane
  Y = reshape(Y(3:5:end,:)/M0_,[],numel(tspan));
  plot(Y(1,:),Y(2,:),'k','LineWidth',1);
end

% critical points
plot(x0(i),x0(i),'r^','MarkerFaceColor','r','MarkerSize',10);
plot(x1(1,i),x2(1,i),'bo','MarkerFaceColor','b','MarkerSize',10);
plot(x1(2,i),x2(2,i),'b^','MarkerFaceColor','b','MarkerSize',10);

% make this plot publishable
set(gca,'TickLabelInterpreter','latex');
set(gca,'xtick',[0 1],'ytick',[1],'fontsize',9);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[200 200 260 260]);
%print -depsc ~/Desktop/3cell_phase_plot.eps

return;

%% (2) 2-cell problem
W = [0 1; 0.5 0.5];

% sweep
lam = nan(3,numel(svec));
res = nan(3,numel(svec));
for i = 1:numel(svec)
  s = svec(i);

  % reduced parameters (a,b)
  [a,b] = hes1reduce(alpha./[1 s 1 1 1],mu,H);

  % useful defines
  f = @(x)1./(a+x.^k);
  g = @(x)1./(1+b*x.^h);
  df = @(x)-k*x.^(k-1)./(a+x.^k).^2;
  dg = @(x)-b*h*x.^(h-1)./(1+b*x.^h).^2;
  F = @(x)(W*g(x)).*f(x)-x;
  J = @(x)diag(f(x))*W*diag(dg(x))+diag(W*g(x).*df(x))-eye(2);

  % spectra
  xroot = [x0(i) x0(i)]';
  J_ = J(xroot);
  lam(1,i) = max(real(eig(J_)));
  res(1,i) = norm(F(xroot));
  xroot = [x1(1,i) x2(1,i)]';
  J_ = J(xroot);
  lam(2,i) = max(real(eig(J_)));
  res(2,i) = norm(F(xroot));
  xroot = [x1(2,i) x2(2,i)]';
  J_ = J(xroot);
  lam(3,i) = max(real(eig(J_)));
  res(3,i) = norm(F(xroot));
end

figure(3), clf, hold on,
ix = find(lam(2,:) > 0);
plot(svec(ix),[x1(1,ix); x2(1,ix)],'b--', ...
     'LineWidth',2,'HandleVisibility','off');
ix = find(lam(2,:) <= 0);
plot(svec(ix),x1(1,ix),'b', ...
     'LineWidth',2,'HandleVisibility','off');
plot(svec(ix),x2(1,ix),'b', ...
     'LineWidth',2);

ix = find(lam(3,:) > 0);
plot(svec(ix),[x1(2,ix); x2(2,ix)],'b--', ...
     'LineWidth',2,'HandleVisibility','off');
ix = find(lam(3,:) <= 0);
plot(svec(ix),[x1(2,ix); x2(2,ix)],'b', ...
     'LineWidth',2,'HandleVisibility','off');

ix = find(lam(1,:) > 0);
plot(svec(ix),x0(ix),'r--','LineWidth',2, ...
     'HandleVisibility','off')
ix = find(lam(1,:) <= 0);
plot(svec(ix),x0(ix),'r','LineWidth',2);
xline(scrit,'k','HandleVisibility','off');
xline(scrit2,'k-.','HandleVisibility','off');

%% 2-cell phase portrait
i = 1;
s = svec(i);

% reduced parameters (a,b)
[a,b] = hes1reduce(alpha./[1 s 1 1 1],mu,H);
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b*x.^h);
  
X1 = linspace(0,1,41);
X2 = linspace(0,1,39)';
[X1_,X2_] = meshgrid(X1,X2);
Y1 = g(X1);
Y2 = g(X2);
F1 = Y2.*f(X1)-X1;
F2 = (0.5*Y1+0.5*Y2).*f(X2)-X2;
D = sqrt(F1.^2+F2.^2);

figure(4), clf,
quiver(X1_,X2_,F1,F2,5);
hold on,
plot(x0(i),x0(i),'ro','MarkerFaceColor','r');
plot(x1(1,i),x2(1,i),'bo','MarkerFaceColor','b');
plot(x1(2,i),x2(2,i),'ko','MarkerFaceColor','k');
contour(X1_,X2_,D,10.^-[0 0.25 0.5 1 2 3 4 5 6 7 8],'b')
axis([0 1 0 1]);

% ----------------------------------------------------------------------
function res = l_criterion(alpha,mu,H)
%L_CRITERION Criterion for stability of homogeneous solution.

[a,b] = hes1reduce(alpha,mu,H);
k = H(3); h = H(4);
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b.*x.^h);
df = @(x)-k*x.^(k-1).*f(x).^2;
dg = @(x)-h*b.*x.^(h-1).*g(x).^2;
x0 = hes1red_homogeneous(a,b,k,h);
res = f(x0).*dg(x0)*0.5-df(x0).*g(x0)+1;

end
% ----------------------------------------------------------------------
function res = l_criterion2(alpha,mu,H)
%L_CRITERION2 Criterion for existence of non-homogeneous solutions.

[a,b] = hes1reduce(alpha,mu,H);
k = H(3); h = H(4);
f = @(x)1./(a+x.^k);
g = @(x)1./(1+b.*x.^h);
df = @(x)-k*x.^(k-1).*f(x).^2;
dg = @(x)-h*b.*x.^(h-1).*g(x).^2;

[~,x12] = hes1red_nonhomogeneous(a,b,k,h);
if any(isnan(x12))
  res = 1;
  return;
end
x1 = x12(1); x2 = x12(2);
res = 0.5*dg(x1)*f(x2)/ ...
      (1-0.5*(df(x2)*g(x2)+f(x2)*dg(x2))-0.5*g(x1)*df(x2))* ...
      f(x1)*dg(x2)/(1-df(x1)*g(x2))-1;

end
% ----------------------------------------------------------------------
