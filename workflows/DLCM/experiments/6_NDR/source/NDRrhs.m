function dX = NDRrhs(t,X,Nj,Np,param)
%NDRrhs Right-hand side for Notch-Delta-Reporter model.

% S. Engblom 2018-02-08 (Minor revision)
% S. Engblom 2017-09-11 (Revision)
% S. Engblom 2017-01-17
% S. Engblom 2016-11-24

% X = [N; D; R], vectorized syntax
len = size(X,1)/3;
N = X(1:len,:);
D = X(len+1:2*len,:);
R = X(2*len+1:end,:);

% (1) sum of junctional contacts
N_a = Nj*N-N; % (note: diagonal convention)
D_a = Nj*D-D;

% (2) sum of protrusional contacts
N_b = Np*N-N;
D_b = Np*D-D;

% (3) form weighted signals
N_in = param.wa*N_a+param.wb*N_b;
D_in = param.wa*D_a+param.wb*D_b;
D_out = param.qa*D_a+param.qb*D_b;

% The Delta-Notch-Reporter model is
%   N' = betaN-<D_in>*N/kt-D*N/kc-gammaN*N
%   D' = betaD*G(R)-D*<N_in>/kt-D*N/kc-gammaD*D
%   R' = betaR*F(<D_out>*N)-gammaR*R
% with gammaN = gammaD = gammaR = 1, and in terms of
%   G(x) = 1/(1+x^m), (m = 2),
%   F(x) = x^s/(kRS+x^s), (s = 2).
%
% The average signals <.> are given by
%   <N_in>  = wa*<N_a>+wb*<N_b>
%   <D_in>  = wa*<D_a>+wb*<D_b>
%   <D_out> = qa*<D_a>+qb*<D_b>
% where <._a> denote junctional contact average and <._b> protrusional
% contact average.

prod = D.*N/param.kc;
dN = param.betaN-D_in.*N/param.kt-prod-N;
dD = param.betaD./(1+R.^2)-D.*N_in/param.kt-prod-D;
DN2 = (D_out.*N).^2;
dR = param.betaR*DN2./(param.kRS+DN2)-R;
dX = [dN; dD; dR];
