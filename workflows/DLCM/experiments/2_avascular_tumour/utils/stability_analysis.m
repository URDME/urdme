%STABILITY_ANALYSIS Stability analysis for the avascular tumor model.
% Sections:
%    1 - Dirichlet BC for pressure at tumor boundary
%    2 - Velocity compatability at tumor-external tissue interface
%    3 - Oxygen analysis (induced perturbations & 'radial' stability analysis)

% E. Blom 2022-11-08 (combined all stability scripts into one)
% S. Engblom 2022-05-30

% *** edit here: what sections to run
sections = [1 2 3 4]

if any(sections == 1)
  disp('*** Section 1 ***');
  
%% Pressure & Dispersion Relation I
% Dirichlet BC for p at tumor boundary, p0 - constant
syms rho_p rho_d rho_pD rho_dD D p0 sigma C P0
syms r r_p r_q r_n
syms c1 c2 d1 d2 e1 e2

% general solutions with different source terms
% Pout = 0
Pp = -rho_pD*r^2/4+d1*log(r)+d2;
Pq = d1*log(r)+d2;
Pn = +rho_dD*r^2/4+e1*log(r)+e2;
% (rho_pD = rho_p/D, rho_dD = rho_d/D)

% solve
e1 = 0; % regularity

% gradient continuity at r_n and r_q
d1 = eval(solve(subs(diff(Pn,r)-diff(Pq,r),r,r_n),d1));
d1 = eval(solve(subs(diff(Pq,r)-diff(Pp,r),r,r_q),d1));

% outer boundary condition = P0 = p0-sigma*C
d2 = eval(solve(subs(Pp,r,r_p)-P0,d2));

% continuity at r_q and r_n
d2 = eval(solve(subs(Pq-Pp,r,r_q),d2));
e2 = eval(solve(subs(Pn-Pq,r,r_n),e2));

% "answer":
%pretty(collect(eval(Pp),r))

% perturbation ansatz
% ~P = P(r) + eps * Q *exp(lambda*t)*cos(k*theta)
% ~rq = rq + eps * rqc *exp(lambda*t)*cos(k*theta) 
% ~rn = rn + eps * rnc *exp(lambda*t)*cos(k*theta) 
% rqc, rnc implied from how oxygen is perturbed by ~rp
syms k C1 C2 D1 D2 E1 E2 rqc rnc
% Qout = 0
Qp = C1*r^k+C2*r^(-k); 
Qq = D1*r^k+D2*r^(-k);
Qn = E1*r^k+E2*r^(-k);

% solve
E2 = 0; % regularity

% 1st order continuity at r_n and r_q
cont = {eval(subs(Qn-Qq,r,r_n)) ...
        eval(subs(Qq-Qp,r,r_q)) ...
        eval(subs(diff(Qn,r)-diff(Qq,r) + ...   % derivative jump
        rnc*(diff(Pn,2,r) - diff(Pq,2,r)),r,r_n)) ...
        eval(subs(diff(Qq,r)-diff(Qp,r) + ...   % derivative jump
        rqc*(diff(Pq,2,r) - diff(Pp,2,r)),r,r_q))};

% final surface tension boundary condition
sol = solve(cont{:},subs(Qp-sigma*(k^2-1)/r^2+eval(diff(Pp,r)),r,r_p), ...
            C1,C2,D1,D2,E1);
C1 = sol.C1;
C2 = sol.C2;
D1 = sol.D1;
D2 = sol.D2;
E1 = sol.E1;

% "answer":
%pretty(eval(Qp))

% dispersion relation
syms LAM

LAM = -D*eval(subs(diff(Pp,2,r)+diff(Qp,r),r,r_p));

% "answer":
pretty(LAM)
%pretty(subs(LAM,k,1))
%pretty(subs(LAM,k,2))

end

if any(sections == 2)
  disp('*** Section 2 ***');

%% Pressure & Dispersion Relation II
% Pressure propagates in the external tissue with Darcy constant D_out
syms rho_p rho_d rho_pD rho_dD D D_out p_out sigma C P0
syms r R r_p r_q r_n    
syms b1 b2 c1 c2 d1 d2 e1 e2
% R is 'far-away' where we set arbitrary P0 for Pout
% Both R and P0 are placeholder variables for below symbolic evaluation
% - they are not used in the final expression!

% far-field pressure analysis (comment away for general condition)
%D = D_out;

% general solutions with different source terms
Pout = b1*log(r) + b2;
Pp = -rho_pD*r^2/4+d1*log(r)+d2;
Pq = d1*log(r)+d2;
Pn = +rho_dD*r^2/4+e1*log(r)+e2;
% (rho_pD = rho_p/D, rho_dD = rho_d/D)

% solve
e1 = 0; % regularity

% gradient continuity at r_n and r_q
d1 = eval(solve(subs(diff(Pn,r)-diff(Pq,r),r,r_n),d1));
d1 = eval(solve(subs(diff(Pq,r)-diff(Pp,r),r,r_q),d1));
b1 = eval(solve(subs(D*diff(Pp,r)-D_out*diff(Pout,r),r,r_p),b1));

% outermost BC: Pout = P0
b2 = eval(solve(subs(Pout,r,R)-P0,b2));

% tumor-external tissue boundary condition: Pp = Pout-sigma*C
d2 = eval(solve(subs(Pp,r,r_p)-(subs(Pout,r,r_p)-sigma*C),d2));

% continuity at r_q and r_n
d2 = eval(solve(subs(Pq-Pp,r,r_q),d2));
e2 = eval(solve(subs(Pn-Pq,r,r_n),e2));

% "answer":
% pretty(collect(eval(Pp),r))

% perturbation ansatz
% ~P = P(r) + eps * Q *exp(lambda*t)*cos(k*theta)
% ~rq = rq + eps * rqc *exp(lambda*t)*cos(k*theta) 
% ~rn = rn + eps * rnc *exp(lambda*t)*cos(k*theta) 
% rqc, rnc implied from how oxygen is perturbed by ~rp
syms k B1 B2 C1 C2 D1 D2 E1 E2 rnc rqc rpc
Qout = B1*r^k+B2*r^(-k);
Qp = C1*r^k+C2*r^(-k);
Qq = D1*r^k+D2*r^(-k);
Qn = E1*r^k+E2*r^(-k);

% solve
E2 = 0; % regularity

% perturbation finite far away
B1 = 0;

% 1st order continuity at r_n and r_q
cont = {eval(subs(Qn-Qq,r,r_n)) ...
        eval(subs(Qq-Qp,r,r_q)) ...
        eval(subs(diff(Qn,r)-diff(Qq,r) + ...           % derivative jump
        rnc*(diff(Pn,2,r) - diff(Pq,2,r)),r,r_n)) ...
        eval(subs(diff(Qq,r)-diff(Qp,r) + ...    % derivative jump
        rqc*(diff(Pq,2,r) - diff(Pp,2,r)),r,r_q)) ...
        eval(subs(D*(diff(Qp,r) + rpc*diff(Pp,2,r)) - ...
        D_out*(diff(Qout,r) + rpc*diff(Pout,2,r)),r,r_p))};

% final surface tension boundary condition
sol = solve(cont{:},subs(Qp-Qout-rpc*sigma*(k^2-1)/r^2+rpc*eval(diff(Pp,r))- ...
                  rpc*eval(diff(Pout,r)),r,r_p),B2,C1,C2,D1,D2,E1);

B2 = sol.B2;
C1 = sol.C1;
C2 = sol.C2;
D1 = sol.D1;
D2 = sol.D2;
E1 = sol.E1;
% rnc, rqc are solved for in Section 3
% rpc is intentionally arbitrary

% "answer":
% pretty(eval(Qp))

% dispersion relation
syms LAM

LAM = -D*eval(subs(rpc*diff(Pp,2,r)+diff(Qp,r),r,r_p))/rpc;

% "answer":
pretty(LAM)
% pretty(subs(LAM,k,1))
% pretty(subs(LAM,k,2))

end

if any(sections == 3)
  disp('*** Section 3 ***');
  
%% Oxygen Perturbation & Regional Sizes w.r.t Time
% Dirichlet BC at R for oxygen. 
% i) Evaluates how boundary perturbation
% propagates to the oxygen regions 

% E. Blom 2022-06-14

syms rho_p rho_d kappa_p kappa_d cons cons_prol c0 
syms r R r_p r_q r_n
syms a1 a2 b1 b2 c1 c2 d1 d2
% cons_prol = cons*(1+delta), delta between 0 and 1 reasonably

% general solutions with different source terms
Cout = a1*log(r)+a2;
Cp = cons_prol*r^2/4+b1*log(r)+b2;
Cq = cons*r^2/4+c1*log(r)+c2;
Cn = d1*log(r)+d2;

% solve
d1 = 0; % regularity
%c1 = 0;    % for rn = 0

% gradient continuity at r_n, r_q, and r_p
c1 = eval(solve(subs(diff(Cn,r)-diff(Cq,r),r,r_n),c1));
b1 = eval(solve(subs(diff(Cq,r)-diff(Cp,r),r,r_q),b1));
a1 = eval(solve(subs(diff(Cp,r)-diff(Cout,r),r,r_p),a1));

% outer boundary condition = C(R) = c0
a2 = eval(solve(subs(Cout,r,R)-c0,a2));

% continuity at r_p, r_q, and r_n
b2 = eval(solve(subs(Cp-Cout,r,r_p),b2));
c2 = eval(solve(subs(Cq-Cp,r,r_q),c2));
d2 = eval(solve(subs(Cn-Cq,r,r_n),d2));

% "answer":
% pretty(collect(eval(cpq),r))

% Regional sizes defined by oxygen thresholds
kappa_p = subs(Cq,r,r_q);
kappa_d = subs(Cq,r,r_n);

% Induced perturbation analysis - finding rnc and rqc

% perturbation ansatz
% ~C = C(r) + eps * Q *exp(lambda*t)*cos(k*theta)
% ~rq = rq + eps * rqc *exp(lambda*t)*cos(k*theta) 
% ~rn = rn + eps * rnc *exp(lambda*t)*cos(k*theta) 
% rqc, rnc implied from how oxygen is perturbed by ~rp
syms k B1 B2 C1 C2 D1 D2 E1 E2 rnc rqc rpc
Qout = B1*r^k+B2*r^(-k);
Qp = C1*r^k+C2*r^(-k);
Qq = D1*r^k+D2*r^(-k);
Qn = E1*r^k+E2*r^(-k);

% solve
E2 = 0; % regularity
%D2 = 0;

% 1st order continuity at r_n and r_q
cont = {eval(subs(Qn-Qq,r,r_n)) ...
        eval(subs(Qq-Qp,r,r_q)) ...
        eval(subs(Qp-Qout,r,r_p)) ...
        eval(subs(diff(Qn,r)-diff(Qq,r) + ...       % derivative jump
        rnc*(diff(Cn,2,r) - diff(Cq,2,r)),r,r_n)) ...
        eval(subs(diff(Qq,r)-diff(Qp,r) + ...  % derivative jump (if delta ~= 0)
        rqc*(diff(Cq,2,r) - diff(Cp,2,r)),r,r_q)) ...  
        eval(subs(diff(Qp,r)-diff(Qout,r) + ...     % derivative jump
        rpc*(diff(Cp,2,r) - diff(Cout,2,r)),r,r_p)) ...
        eval(subs(rnc*diff(Cq,r) + Qq, r, r_n)) ...  % oxygen thresholds
        eval(subs(rqc*diff(Cq,r) + Qq, r, r_q))};

% final surface tension boundary condition
sol = solve(cont{:}, subs(Qout,r,R) ...       % R unperturbed
            ,B1,B2,C1,C2,D1,D2,E1, rqc, rnc);

B1 = sol.B1;
B2 = sol.B2;
C1 = sol.C1;
C2 = sol.C2;
D1 = sol.D1;
D2 = sol.D2;
E1 = sol.E1;
rnc = sol.rnc;
rqc = sol.rqc;

%pretty(rnc)
pretty(rqc)

% Simple tests:
% if cons_prol = cons
% should have pretty(eval(Cq-Cp)) = 0 

% Radial stability analysis INVALID
% 
% % Differentiate (implicitly) these constant kappa:s with respect to time
% % (i.e., 0 = null1 = null2 = null3)
% syms rp(t) rq(t) rn(t)  % introduce t-dependent variables
% 
% %null1 = diff(subs(eval(kappa_p),[r_p, r_q, r_n], [rp, rq, rn]),t);
% % null1 is complicated expression, below two are simpler & equivalent
% 
% null2 = diff(subs(eval(kappa_d),[r_p, r_q, r_n], [rp, rq, rn]),t);
% null3 = diff(subs(eval(kappa_p-kappa_d),[r_p, r_q, r_n], [rp, rq, rn]),t);
% 
% % uncomment to visualise
% %pretty(null2)
% %pretty(null3)

end