%RAD_STABILITY Stability of reduced tumor model.
%   This script investigates the stability conditions for a radially
%   symmetric avascular tumor model. Additionally, the script contains
%   some analysis of inversion formulas related to the reduced tumor
%   model.

% S. Engblom 2023-05-15 (Revision, partially new notation)
% S. Engblom 2023-04-06

% Firstly,
%   K_prol := 4*(1-kappa_prol)/lambda >= 0
% since c = 1 at the oxygen boundary,
%   K_death:= 4*(1-kappa_death)/lambda >= 0
% for the same reason, and also K_prol < K_death since kappa_death <
% kappa_prol.
%
% To sum up 0 < K_prol < K_death by non-dimensionalization and
% physical reasons.
%
% Below we find the additional natural conditions that K_death <= 1
% and exp(-1) <= K_prol <= 1, i.e., all in all exp(-1) <= K_prol <
% K_death <= 1.

%% Let D(r2) = r2*log(r2)-r2. We note that D(r2) < 0 and D'(r2) < 0
% for all 0 < r^2 =: r2 < 1.
%
% The equation C+D(r2) = 0 can be solved for r2 as follows:
%   C+r2*log(r2)-r2 = 0,
% or, after excluding r2 = 0,
%   C/r2+log(r2)-1 = 0.
% Changing sign and taking exponentials we get
%   exp(-C/r2)*1/r2 = exp(-1),
% or
%   exp(-C/r2)*(-C/r2) = -C*exp(-1),
% which has the real solutions
%   -C/r2 = lambertw({0,-1},-C*exp(-1)),
% or using x/lambertw(x) = exp(lambertw(x)),
%   r2 = exp(lambertw({0,-1},-C*exp(-1))+1),
% provided 0 <= C <= 1 (or else there are only complex solutions). The
% branch k = -1 yields a smaller solution for r2 than the k = 0
% branch.

%% The first equation to solve (for rn2) is
%   -K_death-D(rp2)+D(rn2) = 0,
% where 0 < rn2 < rp2. Attempting rn2 = rp2 we get -K_death < 0, while
% for rn2 = 0 we get the requirement that
%   D2 := -K_death-D0 := -K_death-D(rp2) >= 0, (1)
% or else there is no solution (by continuity and monotonicity).
%
% The solutions is then
%   rn2 = exp(lambertw(-1,-D2*exp(-1))+1),
% This adds the condition that
%   D2 <= 1, (2)
% or else the solution is complex-valued.

% To conclude, rn2 can be solved according to the above, with the
% conditions on D2 given as (1) and (2) above or else there are no
% solutions.
%
% The necrotic phase starts for (by solving (1))
%   rp2 = exp(lambertw(-1,-K_death*exp(-1))+1). [A]
% In particular, if K_death > 1, then there are no solutions for any
% value of rp (no necrotic region). We also have for K_death = 0 that
% rp^2 = 0 such that the whole tumor becomes necrotic in this limit.

%% Let E(r2) = B*log(r2)-r2 for 0 <= B < r^2 (actually, B = rn^2)
% such that again E(r2) < 0 and E'(r2) < 0 for the r under
% consideration.
%
% The equation C+E(r2) = 0 is solved as follows:
%   C+B*log(r2)-r2 = 0,
%   exp(C)*r2^B*exp(-r2) = 1,
%   (-r2/B)*exp(-r2/B) = -1/B*exp(-C/B),
% which has the real solutions
%  -r2/B = lambertw({0,-1},-1/B*exp(-C/B))
% provided 0 > -1/B*exp(-C/B) >= -exp(-1). Here the branch k = -1
% yields the larger solution for r2 and hence r = r^2.

%% The second equation to solve (for rq2 now) is
%   -K_prol-D(rp2)+E(rq2) = 0,
% where E(rq2) = rn2*log(rq2)-rq2 and where 0 < rn2 < rq2 < rp2. For
% rq2 = rn2 we find -K_prol-D(rp2)+D(rn2) > -K_death-D(rp2)+D(rn2) = 0,
% while for rq2 = rp2 we have
%   D1 := -K_prol-D0 := -K_prol-D(rp2) < -E(rp2) < 1/2, (3)
% or else there is no solution (by continuity and monotonicity).
%
% The solution is thus
%   rq2 = -rn2*lambertw(-1,-1/rn2*exp(-D1/rn2)).
% Note that for rn2 = 0, rq2(rp2) = D1 trivially. In this case we find
% that the quiescent phase starts when
%   rp2 = exp(lambertw(-1,-K_prol*exp(-1))+1). [B]
% The formula is identical to [A] but with K_death replaced with
% K_prol and so the same conclusions as before hold. In particular, if
% K_prol > 1, then there will not be a quiescent region.
%
% In fact, if
%   D1 = rp2
% can be solved wrt to rp3, then this implies a stationary value
% without a necrotic region. This occurs for K_prol <= exp(-1) and
% yields
%   rp2 = -K_prol/lambertw(-1,-K_prol) = exp(lambertw(-1,-K_prol)).

%% Test all of this numerically.

if ~exist('ncase','var'), ncase = 1; end 
if ~exist('noplot','var'), noplot = false; end 

if ncase ~= 0 && ncase ~= -1
  % select 0 < K_prol < K_death <= 1 first
  PAR = {[0.40 0.80 5] ...
         [0.36 0.80 5] ...    % K_prol < exp(-1), stationary before necrotic
         [0.40 0.80 4.75] ... % no stationary point
        };
  K_prol = PAR{ncase}(1); % 4*(1-kappa_prol)/lambda
  K_death = PAR{ncase}(2); % 4*(1-kappa_death)/lambda
  mu_death = PAR{ncase}(3);
elseif ncase == -1
  [K_prol,K_death,mu_death,discr] = l_inveq(rn_eq^2,rq_eq^2,rp_eq^2);
end

%D = @(r2)r2.*log(r2)-r2;

% left limits for the necrotic and quiescent phases (<= 1 trivially)
rp_min_n = real(sqrt(exp(lambertw(-1,-K_death*exp(-1))+1)));
rp_min_q = real(sqrt(exp(lambertw(-1,-K_prol*exp(-1))+1)));

rp = linspace(0,1,1e3);
rp2 = rp.^2;
[rn2,rq2] = l_rads2(rp2,K_prol,K_death);
rn = sqrt(rn2);
rq = sqrt(rq2);

% find equilibrium rp = rp_max (if any)
rp_max = find(mu_death*rn2+rq2-rp2 >= 0,2,'first');
if ~isempty(rp_max), rp_max = rp_max(2:end); end % (avoid rp = 0)
if ~isempty(rp_max)
  % solve
  rp2_max = fsolve(@(rp2)l_Feq(rp2,K_prol,K_death,mu_death), ...
                  rp2(rp_max(1)),optimoptions('fsolve','Display','none'));
  rp2_max = real(rp2_max);
  rp_max = sqrt(rp2_max);
end

if ~noplot
  figure(1), clf,
  plot(rp,rp,'r.-', ...
       rp,rn,'k.-', ...
       rp,rq,'g.-');
  hold on,
  xline(rp_min_n,'k');
  xline(rp_min_q,'g');
  if ~isempty(rp_max), xline(rp_max,'r--'); end
  legend('r_p','r_n','r_q','location','nw');
  xlabel('r_p');
  switch ncase,
    case  1, title('All regions exist');
    case  2, title('Necrotic region does not exist');
    case  3, title('No stationary point');
    case  0, title(sprintf('C1 = %f, C2 = %f, mu\\_death = %f',C1,C2,mu_death));
    case -1, title(sprintf(['rn\\_eq = %f, rq\\_eq = %f, rp\\_eq = ' ...
                          '%f, mu\\_death =  = %f'], ...
                           rn_eq,rq_eq,rp_eq,mu_death));
  end
end

if ~isempty(rp_max), 
  [rn2_max,rq2_max] = l_rads2(rp2_max,K_prol,K_death);
  discr = l_discriminant(rn2_max,rq2_max,rp2_max,mu_death);
end
if ~noplot
  figure(2), clf,
  sign_sqrt = @(x)sign(x).*sqrt(abs(x));
  plot(rp,sign_sqrt(l_Feq(rp2,K_prol,K_death,mu_death)),'linewidth',2);
  hold on,
  plot(rp,l_discriminant(rn2,rq2,rp2,mu_death),'linewidth',2);
  xlabel('r_p');
  title('RHS and Discriminant');
  xline(rp_min_n,'k');
  xline(rp_min_q,'g');
  if ~isempty(rp_max), 
    xline(rp_max,'r--');
    plot(rp_max,0,'ro');
    plot(rp_max,discr,'r*');
  end
  yline(0,'k--'), grid on
  
  figure(3), clf,
  [k_prol,k_death] = l_invariant(rn2,rq2,rp2);
  plot(rp,k_prol,rp,k_death,'linewidth',2);
  yline([K_prol,K_death],'k--','linewidth',2);
  title('Invariants');
  xline(rp_min_n,'k');
  xline(rp_min_q,'g');
  legend('K_{prol}','K_{death}');
  xlabel('r_p');
end

return;

%% Proposition 3.1 Tests

% discriminant when rp_eq and mu_d are given
theta = linspace(0,1,1000); theta([1 end]) = [];
% assumed << 1, but > exp(-1) ~ 0.368 or else always stable
rp_eq =  0.4;
mu_d = 1e-5; % assumed > 0
rq_eq = rp_eq./sqrt(1+mu_d*theta);
rn_eq = sqrt(theta).*rp_eq;
discr = l_discriminant(rn_eq.^2,rq_eq.^2,repmat(rp_eq.^2,size(theta)),mu_d);

figure(4), clf,
plot(theta,discr); hold on,
yline(0);
xlabel('\theta'); ylabel('discr');

% the monotonicity of lambda_r w.r.t. theta, while mudeath = 0
rp = [0.6 0.5 0.4 exp(-1)]';
theta = linspace(0,1,1000);
discr = 1+log(rp)./(log(rp)+0.5*log(theta)).*(log(theta)./(1-theta));

figure(5), clf,
plot(theta,discr,'linewidth',2);
hold on,
yline(0,'r--','linewidth',2)
legend(num2str(rp),'location','SE');
ylabel('$\lambda_r$', 'Interpreter','latex')
xlabel('$\theta$', 'Interpreter','latex')
title('$\lambda_r(\theta, \mu_{death} \to 0+)$ for different $r_p^{eq}$',...
    'Interpreter','latex')

% lambda_r for rp = exp(-1) ~ 0.368
figure(6), clf,
discr = 1+1./(1-0.5*log(theta)).*log(theta)./(1-theta);
plot(theta,discr,'linewidth',2);
ylabel('$\lambda_r$', 'Interpreter','latex')
xlabel('$\theta$', 'Interpreter','latex')

% elementary inequality
figure(7), clf,
y = linspace(0,10,1000);
plot(y,-y,'b',y,(exp(-y)-1).*(1+y/2),'r');

return;

%% a "game" in (rn,rq,rp)_eq-space, copy and paste this and then enjoy
ncase = -1;
PAR = {linspace(0,1,101) ...
       linspace(0,1,101) ...
       linspace(0,1,101)};
ii = 6;
jj = 15;
kk = 30;
while 1
  rn_eq = PAR{1}(ii);
  rq_eq = PAR{2}(jj);
  rp_eq = PAR{3}(kk);
  rad_stability;
  key = input('type 1..6 or 0 for exit, then hit enter > ');
  switch key
    case 1, ii = max(ii-1,1);
    case 2, ii = min(ii+1,numel(PAR{1}));
    case 3, jj = max(jj-1,1);
    case 4, jj = min(jj+1,numel(PAR{2}));
    case 5, kk = max(kk-1,1);
    case 6, kk = min(kk+1,numel(PAR{3}));
    case 0, disp('Thanks for playing!'); break;
    otherwise, warning('Only 1..6 or 0!');
  end
end

% ----------------------------------------------------------------------
function F = l_Feq(rp2,K_prol,K_death,mu_death)
%L_FEQ RHS equation of equilibrium.

[rn2,rq2] = l_rads2(rp2,K_prol,K_death);
F = mu_death*rn2+rq2-rp2;

end
% ----------------------------------------------------------------------
function D = l_discriminant(rn2,rq2,rp2,mu_death)
%L_DISCRIMINANT Discriminant.
%   Note: assumes all r-vectors to have the same sizes.

ix0 = abs(rn2) <= 1e-8; % rn ~ 0
ix1 = ~ix0 & abs(rn2-rq2) <= 1e-8*rq2; % rn ~ rq
ix = ~ix0 & ~ix1;
D = zeros(size(rp2));
D(ix0) = 1+log(rp2(ix0));
D(ix1) = 1-log(rp2(ix1))./log(rn2(ix1)).*(mu_death+1);
D(ix) = 1-log(rp2(ix))./log(rn2(ix)).* ...
        (mu_death+log(rq2(ix)./rn2(ix))./(rq2(ix)-rn2(ix)).*rq2(ix));

end
% ----------------------------------------------------------------------
function [K_prol,K_death,mu_death,discr] = l_inveq(rn2_eq,rq2_eq,rp2_eq)
%L_INVEQ Determine parameters for and stability of equilibrium.

D = @(r2)r2.*log(r2)-r2;
D0 = D(rp2_eq);

% immediate: invariants
K_death = -D0+D(rn2_eq);
if rn2_eq == 0, K_death = -D0; end
K_prol = -D0+rn2_eq.*log(rq2_eq)-rq2_eq;
if rn2_eq == 0, K_prol = -D0; end

% immediate: mu_death & discriminant
mu_death = (rp2_eq-rq2_eq)/rn2_eq;
discr = l_discriminant(rn2_eq,rq2_eq,rp2_eq,mu_death);

end
% ----------------------------------------------------------------------
function [k_prol,k_death] = l_invariant(rn2,rq2,rp2)
%L_INVARIANT Invariant relations from square radii.

D = @(r2)r2.*log(r2)-r2;

ix0 = abs(rn2) <= 1e-8; % rn ~ 0
ix = ~ix0;

k_prol(ix0) = -D(rp2(ix0))-rq2(ix0);
k_death(ix0) = -D(rp2(ix0));

k_prol(ix) = -D(rp2(ix))+rn2(ix).*log(rq2(ix))-rq2(ix);
k_death(ix) = -D(rp2(ix))+D(rn2(ix));

end
% ----------------------------------------------------------------------
function [rn2,rq2] = l_rads2(rp2,K_prol,K_death)
%L_RADS2 Determine rn2 and rq2 from rp2 and parameters K_prol, K_death.

D = @(r2)r2.*log(r2)-r2;
D0 = D(rp2);
D2 = -K_death-D0;
D1 = -K_prol-D0;

% rn2 = rn^2 first, avoid complex values
ix = 0 <= D2 & D2 <= 1;
rn2 = zeros(size(rp2)); % ensure 0 <= rn2

% formula:
arg = -D2(ix)*exp(-1);
rn2(ix) = exp(lambertw(-1,arg)+1);
rn2 = real(rn2); % clear roundoff errors

% Identity: exp(W(z)) = z/W(z)
%
% rq2 = -rn2*W(-1/rn2*exp(-D1/rn2))
%
% exp(-D1/rn2) = exp(-D1*exp(-1)*W(X)/X) =
% exp(D1/D2*W(X)) = rn2^(D1/D2)*exp(-D1/D2)
%
% -1/rn2*exp(-D1/rn2) = -rn2^(D1/D2-1)*exp(-D1/D2) = 
% -exp(-1)*Y^(D1/D2-1), Y := rn2*exp(-1)
% -exp(1)*Y*W(-exp(-1)*Y^(D1/D2-1))
%
% = -exp(1)*Y*[-exp(-1)*Y^(D1/D2-1)] * 
%   W(-exp(-1)*Y^(D1/D2-1))/[-exp(-1)*Y^(D1/D2-1)]
% = Y^D1/D2 * exp(-W(-exp(-1)*Y^(D1/D2-1)))
% = exp(1)*Y * exp(-1)*Y^(D1/D2-1) * exp(-W(-exp(-1)*Y^(D1/D2-1)))
% = Y*Z*exp(1-W(-Z)), Z := Y^(D1/D2-1)*exp(-1)

% rq = rq^2, again avoid complex values
ix = rn2 > 0;
rq2 = zeros(size(rp2));
rq2(~ix) = max(D1(~ix),0); % ensure 0 <= rq2

% use expansion for small values of rn2:
ixa = ix & rn2 < D1/100;
ixb = ix & rn2 >= D1/100;

% use the series expansion lambertw(-1,-x) = log(x)-log(-log(x))+R,
% where R ~ log(-log(x))/log(x) ~ <small>...
arg = rn2(ixa).*log(rn2(ixa))+D1(ixa);
rq2(ixa) = D1(ixa)+rn2(ixa).*log(arg);

% second formula:
%rq2(ixb) = -rn2(ixb).*lambertw(-1,-1./rn2(ixb).*exp(-D1(ixb)./rn2(ixb)));
% a bit prettier...
%arg = rn2(ixb)*exp(-1);
%rq2(ixb) = -exp(1)*arg.*lambertw(-1,-exp(-1)*arg.^(D1(ixb)./D2(ixb)-1));
% ...and even somewhat prettier still:
arg1 = rn2(ixb)*exp(-1);
arg2 = arg1.^(D1(ixb)./D2(ixb)-1)*exp(-1);
rq2(ixb) = arg1.*arg2.*exp(1-lambertw(-1,-arg2));
rq2 = real(rq2); % clear roundoff errors

% finalize
rq2 = min(rp2,rq2); % ensure rq <= rp

end
% ----------------------------------------------------------------------
