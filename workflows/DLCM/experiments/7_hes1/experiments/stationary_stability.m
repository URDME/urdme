%STATIONARY_STABILITY Spectral analysis.
%   This script produces spectra for small perturbations from
%   stationary states.
%
%   For homogeneous stationary states, models in, e.g., variables
%   (N1,D1,N2,D2) are rewritten in (N,D,dN,dD) where the latter two
%   are considered small around a perturbation from the homogeneous
%   state. After decoupling and by a linearization, the spectra for
%   this part is derived.
%
%   For non-homogeneous stationary states this trick does not work and
%   the full Jacobian needs to be considered.
%
%   Models tested here are Hes1@2/Hes1@1 and Hes1@5. For the latter,
%   it is checked that the reduction map is positive definite.
%
%   See also STATIONARY_BOUNDS.

% S. Engblom 2024-09-25 (Revision)
% S. Engblom 2024-08-21 (Revision)
% S. Engblom 2024-07-09 (Revision)
% S. Engblom 2024-04-27

%% §1 the Jacobian for small perturbations from
% homogeneous state, Hes1@2 (Proposition 3.2)
syms N1 D1 N2 D2 f(x) g(x) df(x) dg(x) v lam
f1 = D2*f(N1)-N1;
f2 = v*(g(N1)-D1);
f3 = D1*f(N2)-N2;
f4 = v*(g(N2)-D2);

% other patterns:
% pattern = 2:
%f3 = (D1+D2)/2*f(N2)-N2;
% pattern = 4:
%f3 = (D1/3+2*D2/3)*f(N2)-N2;
% pattern = 6:
%f3 = (D1/6+5*D2/6)*f(N2)-N2;

% variable transformation [N D dN dD]' = T*[N1 D1 N2 D2]'
syms N D dN dD
T = [1 0 1 0; 0 1 0 1; 1 0 -1 0; 0 1 0 -1]/2;
F = subs(T*[f1 f2 f3 f4].',[N1 D1 N2 D2].',T\[N D dN dD].');

% keep only up until quadratic in [dN dD]
F = taylor(F,[dN dD],'order',2);

% quick way to get to the Jacobian in [dN dD] only:
dJ = jacobian(F(3:4,:),[dN dD]);
dJ = subs(dJ,[D diff(f(N),N) diff(g(N),N)],[g(N) df(N) dg(N)]);

lambda_Hes1_2 = eig(dJ)

% Hes1@1 - simplifies to scalar
g1 = g(N2)*f(N1)-N1;
g2 = g(N1)*f(N2)-N2;
T = [1 1; 1 -1]/2;
F = subs(T*[g1 g2].',[N1 N2].',T\[N dN].');
F = taylor(F,dN,'order',2);
dJ = jacobian(F(2,:),dN);
dJ = subs(dJ,[D diff(f(N),N) diff(g(N),N)],[g(N) df(N) dg(N)]);
lambda_Hes1_1 = eig(dJ)

%% §2 the Jacobian for small perturbations from
% non-homogeneous state, Hes1@2 (Proposition 3.4)
F = [f1 f2 f3 f4].';
dJ = jacobian(F,[N1 D1 N2 D2]);
dJ = subs(dJ,[D1 diff(f(N1),N1) diff(g(N1),N1)],[g(N1) df(N1) dg(N1)]);
dJ = subs(dJ,[D2 diff(f(N2),N2) diff(g(N2),N2)],[g(N2) df(N2) dg(N2)]);

% cofactor expansion of characteristic polynomial, peel off one factor
% at a time
q = coeffs(charpoly(dJ,lam),lam,'all');

dJ = lam*eye(4)-dJ;
% $$$ dJ(1,1)*det(dJ(2:4,2:4))-dJ(1,4)*det(dJ(2:4,1:3));
% $$$ dJ(1,1)*dJ(2,2)*det(dJ(3:4,3:4))-dJ(1,4)*dJ(2,1)*det(dJ(3:4,2:3));
p = dJ(1,1)*dJ(2,2)*dJ(3,3)*dJ(4,4)-dJ(1,4)*dJ(2,1)*dJ(3,2)*dJ(4,3);

% check that this = 0
q-coeffs(collect(expand(p),lam),lam,'all')

% the 0th order coefficient is the discriminant:
pvec = coeffs(collect(expand(p),lam),lam);
p0 = collect(pvec(1),[v])
% by Descartes' rule of sign there is an eigenvalue with positive real
% part <==> p0 < 0

% Hes1@1
F = [g1 g2].';
dJ = jacobian(F,[N1 N2]);
dJ = subs(dJ,[diff(f(N1),N1) diff(g(N1),N1)],[df(N1) dg(N1)]);
dJ = subs(dJ,[diff(f(N2),N2) diff(g(N2),N2)],[df(N2) dg(N2)]);
p = (lam-dJ(1,1))*(lam-dJ(2,2))-dJ(1,2)*dJ(2,1)
q = coeffs(charpoly(dJ,lam),lam,'all');
q-coeffs(collect(expand(p),lam),lam,'all')
% all coefficients are positive save for the 0th order coefficients
pvec = coeffs(collect(expand(p),lam),lam);
p0 = pvec(1)

%% §3 the Jacobian for small perturbations from
% homogeneous state, Hes1@5
syms f(x) g(x) df(x) dg(x) lam
syms D1 N1 M1 P1 n1 D2 N2 M2 P2 n2
syms aD aN aM aP an mD mN mM mP mn KM Kn
f1 = aD*n1-mD*D1;
f2 = aN*D2-mN*N1;
f3 = aM*N1*f(P1/KM)-mM*M1; % f(x) = 1/(1+x^k)
f4 = aP*M1-mP*P1;
f5 = an*g(P1/Kn)-mn*n1; % g(x) = 1/(1+x^h)
f6 = aD*n2-mD*D2;
f7 = aN*D1-mN*N2;
f8 = aM*N2*f(P2/KM)-mM*M2;
f9 = aP*M2-mP*P2;
f10 = an*g(P2/Kn)-mn*n2;

% variable transformation, e.g., [N dN]' = T*[N1 D1]'
syms D N M P n dD dN dM dP dn 
T = [eye(5) eye(5); eye(5) -eye(5)]/2;
F = subs(T*[f1 f2 f3 f4 f5 f6 f7 f8 f9 f10].', ...
    [D1 N1 M1 P1 n1 D2 N2 M2 P2 n2].',T\[D N M P n dD dN dM dP dn].');

% keep only up until quadratic in small variables dN etc.
F = taylor(F,[dD dN dM dP dn],'order',2);

% check that the variables decouple:
dJ = jacobian(F,[D N M P n dD dN dM dP dn]);
dJ(1:5,6:end) % = 0

% check that the reductions are definite:
det(-dJ(6:8,6:8)) % > 0, Hes1@2 in variables (x,y) = (P,n)
det(-dJ([6:8 10],[6:8 10])) % > 0, Hes1@1 in variable x = P

% quick way to get to the Jacobian in small variables only:
dJ = jacobian(F(6:end,:),[dD dN dM dP dn]);
dJ = subs(dJ,[diff(f(P/KM),P) diff(g(P/Kn),P)],[df(P/KM)/KM dg(P/Kn)/Kn]);

% (note that dJ only depends on the outer variables)

% cofactor expansion of characteristic polynomial, peel off one factor
% at a time
q = coeffs(charpoly(dJ,lam),lam,'all');

dJ = lam*eye(5)-dJ;
% $$$ dJ(1,1)*det(dJ(2:5,2:5))- ...
% $$$     dJ(2,1)*det(dJ([1 3:5],2:5));
% $$$ 
% $$$ dJ(1,1)*dJ(2,2)*det(dJ(3:5,3:5))+ ...
% $$$     dJ(2,1)*dJ(1,5)*det(dJ(3:5,2:4));
% $$$ 
% $$$ dJ(1,1)*dJ(2,2)*dJ(5,5)*det(dJ(3:4,3:4))+ ...
% $$$     dJ(2,1)*dJ(1,5)*dJ(3,2)*det(dJ(4:5,3:4));
p = dJ(1,1)*dJ(2,2)*dJ(5,5)*(dJ(3,3)*dJ(4,4)-dJ(3,4)*dJ(4,3))+ ... 
    dJ(2,1)*dJ(1,5)*dJ(3,2)*dJ(4,3)*dJ(5,4);

% check that this = 0
simplify(q-coeffs(collect(expand(p),lam),lam,'all'))

% the 0th order coefficient is the discriminant:
pvec = coeffs(collect(expand(p),lam),lam);
p0 = collect(pvec(1),[df(P) N aM aP mN mn mM mP])
% by Descartes' rule of sign there is an eigenvalue with positive real
% part <==> p0 < 0

%% §4 the Jacobian for small perturbations from
% non-homogeneous state, Hes1@5
F = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10].';
dJ = jacobian(F,[D1 N1 M1 P1 n1 D2 N2 M2 P2 n2]);
dJ = subs(dJ,[diff(f(P1/KM),P1) diff(g(P1/Kn),P1)],[df(P1/KM)/KM dg(P1/Kn)/Kn]);
dJ = subs(dJ,[diff(f(P2/KM),P2) diff(g(P2/Kn),P2)],[df(P2/KM)/KM dg(P2/Kn)/Kn]);

% check that the reductions are definite:
det(-dJ([1:3 6:8],[1:3 6:8])) % > 0, Hes1@2 in variables (x,y) = (P,n)
det(-dJ([1:3 5 6:8 10],[1:3 5 6:8 10])) % > 0, Hes1@1 in variable x = P

% characteristic polynomial
q = coeffs(charpoly(dJ,lam),lam,'all');

% cofactor expansion of characteristic polynomial, peel off one factor
% at a time
dJ = lam*eye(10)-dJ;
% $$$ p = dJ(1,1)*det(dJ(2:10,2:10))+dJ(7,1)*det(dJ([1:6 8:10],2:10));
% $$$ p = dJ(1,1)*dJ(2,2)*det(dJ(3:10,3:10))- ...
% $$$     aD^2*aM^2*aN^2*aP^2*an^2*dg(P1)*dg(P2)*f(P1)*f(P2);
% $$$ p = dJ(1,1)*dJ(2,2)*dJ(7,7)*det(dJ([3:6 8:10],[3:6 8:10]))- ...
% $$$     aD^2*aM^2*aN^2*aP^2*an^2*dg(P1)*dg(P2)*f(P1)*f(P2);
% $$$ p = dJ(1,1)*dJ(2,2)*dJ(7,7)*dJ(6,6)*det(dJ([3:5 8:10],[3:5 8:10]))- ...
% $$$     aD^2*aM^2*aN^2*aP^2*an^2*dg(P1)*dg(P2)*f(P1)*f(P2);
% $$$ p = dJ(1,1)*dJ(2,2)*dJ(7,7)*dJ(6,6)* ...
% $$$     det(dJ(3:5,3:5))*det(dJ(8:10,8:10))- ...
% $$$     aD^2*aM^2*aN^2*aP^2*an^2*dg(P1)*dg(P2)*f(P1)*f(P2);
% $$$ p = dJ(1,1)*dJ(2,2)*dJ(7,7)*dJ(6,6)*dJ(5,5)*dJ(10,10)* ...
% $$$     det(dJ(3:4,3:4))*det(dJ(8:9,8:9))- ...
% $$$     aD^2*aM^2*aN^2*aP^2*an^2*dg(P1)*dg(P2)*f(P1)*f(P2);
p = dJ(1,1)*dJ(2,2)*dJ(7,7)*dJ(6,6)*dJ(5,5)*dJ(10,10)* ...
    (dJ(3,3)*dJ(4,4)-dJ(4,3)*dJ(3,4))* ...
    (dJ(8,8)*dJ(9,9)-dJ(9,8)*dJ(8,9))- ...
    aD^2*aM^2*aN^2*aP^2*an^2*dg(P1/Kn)*dg(P2/Kn)*f(P1/KM)*f(P2/KM)/Kn^2;

% check that this = 0
simplify(q-coeffs(collect(expand(p),lam),lam,'all'))

% the 0th order coefficient is the discriminant:
pvec = coeffs(collect(expand(p),lam),lam);
p0 = pvec(1)
% by Descartes' rule of sign there is an eigenvalue with positive real
% part <==> p0 < 0

return;

%% §5 simplify the criterion g'f-f'g < -1 but expressed in
% Hes1@5-parameters
syms alpha [5 1]
syms mu [5 1]
syms KM Kn
c1 = alpha(2)*alpha(3)/mu(2);
c2 = (alpha(4)/(mu(4)*KM))^k;
c2 = expand(c2,'IgnoreAnalyticConstraints',true);
c3 = alpha(1)*alpha(5)/mu(5);
c4 = (alpha(4)/(mu(4)*Kn))^h;
c4 = expand(c4,'IgnoreAnalyticConstraints',true);
D0 = c3/mu(1);
M0 = (c1*D0/(mu(3)*c2))^(1/(k+1));
M0 = expand(M0,'IgnoreAnalyticConstraints',true);
a = 1/(c2*M0^k);
a = expand(a,'IgnoreAnalyticConstraints',true);
a = combine(a,'IgnoreAnalyticConstraints',true);
b = c4*M0^h;
b = expand(b,'IgnoreAnalyticConstraints',true);
b = combine(b,'IgnoreAnalyticConstraints',true);

syms P0
x0 = P0*mu(4)/(alpha(4)*M0);
x0 = expand(x0,'IgnoreAnalyticConstraints',true);
x0 = combine(x0,'IgnoreAnalyticConstraints',true);

syms k h
bnd = -b*h*x0^h/(1+b*x0^h)+k*x0^k/(a+x0^k);
bnd1 = -b*h*x0^h/(1+b*x0^h);
bnd1 = expand(bnd1,'IgnoreAnalyticConstraints',true);
bnd1 = combine(bnd1,'IgnoreAnalyticConstraints',true);
bnd1 = simplify(bnd1);

%bnd2 = k*x0^k/(a+x0^k); % does not work
bnd2 = k*x0^(k+1)*(1+b*x0^h); % uses the homogeneous relation

bnd2 = expand(bnd2,'IgnoreAnalyticConstraints',true);
bnd2 = combine(bnd2,'IgnoreAnalyticConstraints',true);
bnd2 = simplify(bnd2);
