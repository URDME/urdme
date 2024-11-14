%PARAMETRIZATION Hes1 model parametrization.

% S. Engblom 2024-09-12

% called from final_parametrization or not
if ~exist('final','var')
  final = false;
end
if ~exist('save_data','var')
  save_data = false;
end

% input: data from perturbed or unperturbed concentrations
if ~final
  DATA = hes1_conc';
else
  NMC = 1000;
  DATA = hes1_conc(NMC);
  alphas = zeros(5,NMC);
end

% input: parameters considered known [/min]
if ~final
  par = hes1_params;
else
  par = hes1_params(NMC);
end
mu = [par.muD' par.muN' par.muM' par.muP' par.mun'];

% input: choice of Hill coefficients
KM = par.KM(1); Kn = par.Kn(1);
k = par.k(1); h = par.h(1);
H = [KM Kn k h];

%% (§1) determine alpha s.t. we agree wrt the homogeneous solution

% solve X122 = DATA wrt log(alpha)
opts = optimset('Display','off','tolx',1e-8,'tolfun',1e-10);

% solve for alphas for data
if ~final
  target = @(log_alpha)mean(l_X122(exp(log_alpha),mu,H),2).'-DATA;
  [log_alpha,~,flag] = fsolve(target,zeros(1,5),opts);
  if flag ~= 1
    disp(sprintf('Convergence issuses, flag = %d.',flag));
  end
  alpha = exp(log_alpha');
else
  log_alpha = zeros(1,5); % used as starting guess
  for i = 1:NMC
    i
    target = @(log_alpha)mean(l_X122(exp(log_alpha),mu(i,:),H),2).'- ...
             DATA(:,i)';
    [log_alpha,~,flag] = fsolve(target,log_alpha,opts);
    if flag ~= 1
      disp(sprintf('Convergence issuses, flag = %d.',flag));
    end
    alphas(:,i) = exp(log_alpha');
  end
end

if final
  if save_data
    save('../data/alpha_parameterization.mat','alphas')
  end
  return;
end

X0 = l_X0(alpha,mu,H);
X122 = l_X122(alpha,mu,H);

% 3-cell coupling
W = ones(3)-eye(3);
W = W./sum(W,2);

% assemble Jacobian function once
[fun,funJ,funC] = hes1_buildODE;

% $$$ % checks:
% $$$
% $$$ % is it really stationary?
% $$$ norm(hes1_System(X122(:),alpha,mu,H,fun,funC,W),inf)
% $$$
% $$$ % did fsolve work wrt DATA?
% $$$ norm(abs((mean(X122,2)'-DATA)./DATA),inf)

% check spectral properties around X0 and X122
J0 = hes1_Jacobian(repmat(X0,1,3),alpha,mu,H,funJ,funC,W);
lam_ = eig(J0);
tau1 = l_escape(lam_)/60 % time until fate decision [h]
per1 = l_period(lam_)/60 % period [h]

%% (2) evaluate response in (alpha,KM,Kn) to improve on spectral properties

if 0
  disp('Computing...');

  % sweep over (KM,Kn)
  KM_ = 10.^linspace(log10(0.25*KM),log10(2*KM),7);
  Kn_ = 10.^linspace(log10(0.25*Kn),log10(2*Kn),6);
  [KM,Kn] = meshgrid(KM_,Kn_);
  tau1_ = zeros(size(KM));
  per1_ = zeros(size(KM));
  alpha_ = zeros([numel(alpha) size(KM)]);
  X0_ = zeros([numel(DATA) size(KM)]);

  ok = 1;
  for i = 1:size(KM,1)
    for j = 1:size(KM,2)
      % initial guess
      if j > 1
        alpha__ = alpha_(:,i,j-1);
      elseif i > 1
        alpha__ = alpha_(:,i-1,j);
      else
        alpha__ = alpha;
      end
      H_ = [KM(i,j) Kn(i,j) H(3:4)];

      target = @(log_alpha)mean(l_X122(exp(log_alpha),mu,H_),2).'-DATA;
      [log_alpha,~,flag] = fsolve(target,log(alpha__'),opts);
      if flag ~= 1
        disp(sprintf('Convergence issuses, flag = %d.',flag));
      end
      alpha__ = exp(log_alpha');
      alpha_(:,i,j) = alpha__;
      X0_(:,i,j) = l_X0(alpha__,mu,H_);
      X122_ = l_X122(alpha__,mu,H_);

      % checks: are we still stationary and ok wrt DATA?
      ok = ok && ...
        norm(hes1_System(repmat(X0_(:,i,j),1,3), ...
        alpha__,mu,H_,fun,funC,W),inf) < 1e-6 && ...
        norm(hes1_System(X122_, ...
        alpha__,mu,H_,fun,funC,W),inf) < 1e-6 && ...
        norm((mean(X122_,2)'-DATA)./DATA,inf) < 1e-6;

      J0_ = hes1_Jacobian(repmat(X0_(:,i,j),1,3),alpha__,mu,H_,funJ,funC,W);
      lam_ = eig(J0_);
      tau1_(i,j) = l_escape(lam_);
      per1_(i,j) = l_period(lam_);
    end
  end
  disp('   ...done.');
  ok

  tau1_ = tau1_/60;
  per1_ = per1_/60;

  figure(1), clf,
  surf(KM,Kn,tau1_);
  hold on,
  surf(KM,Kn,0.5*(6+12)*ones(size(KM))); % target "6--12 hours"
  C = contourc(KM_,Kn_,tau1_-0.5*(6+12)*ones(size(KM)),[0 0]);
  plot3(C(1,2:end),C(2,2:end),0.5*(6+12)*ones(1,size(C,2)-1),'ko-','LineWidth',2);
  plot(C(1,2:end),C(2,2:end),'ko-','LineWidth',2);
  plot(H(1),H(2),'ro','MarkerFaceColor','r');
  xlabel('KM'); ylabel('Kn');
  title('Exit from homogeneous [h]');
  % $$$ [~,ix] = min(abs(tau1_(:)-0.5*(6+12)));
  % $$$ [KM(ix) Kn(ix)] % potentially better values

  figure(2), clf,
  surf(KM,Kn,per1_);
  hold on,
  surf(KM,Kn,0.5*(2+3)*ones(size(KM))); % target "2--3 hours"
  plot3(H(1),H(2),0.5*(2+3),'ro','MarkerFaceColor','r');
  xlabel('KM'); ylabel('Kn');
  title('Period [h]');
  % $$$ [~,ix] = min(abs(per1_(:)-0.5*(2+3)));
  % $$$ [KM(ix) Kn(ix)] % potentially better values

  figure(3), clf,
  surf(KM,Kn,reshape(X0_(4,:,:),size(KM))); hold on,
  plot(H(1),H(2),'ro','MarkerFaceColor','r');
  xlabel('KM'); ylabel('Kn');
  title('P0');

end

% final ODE-solution for the 3-cell problem with chosen parameters
Tend = 24*60*2;
tspan = linspace(0,Tend);

% somewhat carefully pollute the homogeneous solution in the direction
% of the relevant eigenvectors, suitably weighted together
[V,lam_] = eig(J0,'vec');
[~,ix1] = l_escape(lam_);
v1 = mean(V(:,ix1),2); % in case of more than one ix1
v1 = real(v1); v1 = v1/norm(v1);
[~,ix2] = l_period(lam_);
v2 = mean(V(:,ix2),2); % in case of more than one ix2
v2 = real(v2); v2 = v2/norm(v2);
w1 = 0.1; % relative weights
w2 = 1-w1;
weight = 0.1; % total "weight" of pollution
Y0 = repmat(X0,3,1);
Y0 = Y0+weight*norm(Y0)*(w1*v1+w2*v2); assert(all(Y0 > 0));

% solve and visualize
[~,Y] = ode15s(@(t,y)hes1_System(y,alpha,mu,H,fun,funC,W),tspan,Y0, ...
  odeset('AbsTol',1e-8,'RelTol',1e-10));
Y = Y';

figure(4), clf,
semilogy(tspan/60,Y(1:5:end,:),'b','HandleVisibility','off'); hold on,
semilogy(tspan/60,Y(2:5:end,:),'r','HandleVisibility','off');
semilogy(tspan/60,Y(3:5:end,:),'g','HandleVisibility','off');
semilogy(tspan/60,Y(4:5:end,:),'m','HandleVisibility','off');
semilogy(tspan/60,Y(5:5:end,:),'c','HandleVisibility','off');
set(gca,'xtick',0:2:48);
axis([0 tspan(end)/60 1e-5 2]);
xlabel('[h]')

% illustrate some spectral predictions using P as the reference state

% weight*w1*exp(lam1*tau_) = 1, or tau_ = -log(weight*w1)/lam1 =
% -log(weight*w1)*tau1, where tau1 has been computed before
tau0 = -log(weight*w1)*tau1; % (tau1 since we are reasoning around X0)
y0 = max(Y(4:5:end,1:end/2),[],'all');
plot([0 tau0],4*[y0 y0],'k-v');
plot([tau0 tau0],[4*y0 1e-5],'k--','HandleVisibility','off');

% period
[y1,t1] = findpeaks(Y(4,1:end/2),tspan(1:end/2));
plot(t1(1)/60+[0 per1],2*[y1(1) y1(1)],'b-v'); % estimated as per1
plot(t1(1:2)/60,1.5*y1(1:2),'b--.'); % detected via findpeaks
legend(sprintf('Est fate-decision: %2.1fh',tau0), ...
  sprintf('Est period: %2.1fh',per1), ...
  sprintf('Measured period: %2.1fh',(t1(2)-t1(1))/60), ...
  'location','SW');

% ----------------------------------------------------------------------
function [per,ix] = l_period(lam)
%L_PERIOD Estimate period from eigenvalues.
%   PER = L__PERIOD(LAM) computes the period PER from the vector of
%   eigenvalues LAM. The period is defined as the period of the
%   "essentially oscillating" eigenvalue with largest real part. An
%   eigenvalue LAM is deemed essentially oscillating if ABS(IMAG(LAM))
%   >= ABS(REAL(LAM)).
%
%   [PER,IX] = L__PERIOD(LAM) also returns the index into LAM of the
%   corresponding eigenvalue(s). Multiple indices are returned in the
%   case of multiple dominating eigenvalues.

% max (real part) of all essentially oscillating eigenvalues:
ix = find(abs(imag(lam)) >= abs(real(lam)));
if isempty(ix) % case of no essential oscillations:
  per = inf;
  return;
end
lam_ = max(lam(ix),[],'ComparisonMethod','real');

% indices to all with real part close enough:
ix = ix(abs(real(lam(ix))-real(lam_)) < 1e3*eps);
per = 2*pi/mean(abs(imag(lam(ix))));

end
% ----------------------------------------------------------------------
function [tau,ix] = l_escape(lam)
%L_ESCAPE Estimate escape time from eigenvalues.
%   TAU = L_ESCAPE(LAM) estimates the escape time from eigenvalues
%   LAM. The escape time is simply the reciprocal of the largest real
%   part of all eigenvalues.
%
%   [TAU,IX] = L_ESCAPE(LAM) returns indices as above. As with
%   L_PERIOD, multiple indices are returned in the case of multiple
%   dominating eigenvalues.

lam_ = max(lam,[],'ComparisonMethod','real');
if real(lam_) <= 0 % case of no exit:
  tau = inf;
  ix = zeros(0,1);
  return;
end
ix = find((1+1e3*eps)*real(lam) >= real(lam_));
tau = 1/mean(real(lam(ix)));

end
% ----------------------------------------------------------------------
function X0 = l_X0(alpha,mu,H)
%L_X0 Homogeneous steady-state given parameters.
%   Returns X0 = [D0 N0 M0 P0 n0] for given values of alpha_i, mu_i
%   and Hill-functions H = [KM Kn k h], i = [D N M P n].

% reduced scalar model
[a,b,~,~,M0_,D0_] = hes1reduce(alpha,mu,H);

% solve
x0 = hes1red_homogeneous(a,b,H(3),H(4));
y0 = 1/(1+b*x0^H(4));

% translate back into full Hes1@5-space
X0 = hes1unreduce(x0,y0,M0_,D0_,alpha,mu,H,1);

end
% ----------------------------------------------------------------------
function X122 = l_X122(alpha,mu,H)
%L_X122 Non-homogeneous steady-state given parameters.
%   As L_X0 above but returns [X1 X2 X2], the non-homogeneous stable
%   steady-state for the 3-cell problem instead.

% reduced model
[a,b,~,~,M0_,D0_] = hes1reduce(alpha,mu,H);

% assemble full solution
[~,x_] = hes1red_nonhomogeneous(a,b,H(3),H(4));
x0 = [x_ x_(end)]'; % x1 x2 x2
y0 = 1./(1+b*x0.^H(4));

% check
% $$$ F = @(x)[g(x(2))*f(x(1)); 0.5*(g(x(1))+g(x(2)))*f(x(2))]-x'; % = 0
% $$$ F([x0(1) x0(2)])

% translate back into full Hes1@5-space
W = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
X122 = hes1unreduce(x0,y0,M0_,D0_,alpha,mu,H,W);

end
% ----------------------------------------------------------------------
