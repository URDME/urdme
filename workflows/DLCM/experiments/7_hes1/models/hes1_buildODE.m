function [fun,funJ,funC] = hes1_buildODE
%HES1_BUILDODE Hes1 ODE model on a grid.
%   [fun,funJ,funC] = HES1_BUILDODE returns the Hes1 RHS fun, Jacobian
%   funJ, and coupling Jacobian funC as functions. Intended to be
%   called once.
%
%   Examples:
%     [fun,funJ,funC] = hes1_buildODE; % prepare once
%
%     % problem parameters
%     par = hes1_params;
%     W = [0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0]; % 3 cells
%     alpha = [par.alphaD par.alphaN par.alphaM par.alphaP par.alphan];
%     mu = [par.muD par.muN par.muM par.muP par.mun];
%     H = [par.KM par.Kn par.k par.h];
%
%     % (1): sample integration
%     Tend = 24*60*2;
%     tspan = linspace(0,Tend);
%     Y0 = rand(5*size(W,2),1);
%
%     % optional use of Jacobian:
%     opts = odeset('AbsTol',1e-6,'RelTol',1e-8, ...
%        'Jacobian',@(t,y)hes1_Jacobian(y,alpha,mu,H,funJ,funC,W));
%
%     % final solve
%     [~,Y] = ode15s(@(t,y)hes1_System(y,alpha,mu,H, ...
%                            fun,funC,W),tspan,Y0,opts);
%
%     figure(1), clf,
%     plot(tspan/60,Y(:,4:5:end));
%     xlabel('time [min]'); ylabel('P');
%
%     % (2): more advanced integration on a mesh
%     [P,E,T] = basic_mesh(2,12);
%     [V,R] = mesh2dual(P,E,T,'voronoi');
%     [~,~,N] = dt_operators(P,T);
%     N = N./sum(N,2);
%
%     Y0 = rand(5*size(N,2),1);
%     opts = odeset('AbsTol',1e-6,'RelTol',1e-8, ...
%        'Jacobian',@(t,y)hes1_Jacobian(y,alpha,mu,H,funJ,funC,N));
%     [~,Y] = ode15s(@(t,y)hes1_System(y,alpha,mu,H, ...
%                            fun,funC,N),tspan,Y0,opts);
%     Y = reshape(Y',5,[],numel(tspan));
%     Y = squeeze(Y(4,:,:)); % fetch P alone
%
%     % judge fate decision by the end:
%     [Yincr,ix] = sort(Y(:,end));
%     [~,ijmp] = max(diff(Yincr)); % cut at largest jump
%     ixlo = ix(1:ijmp);
%     ixhi = ix(ijmp+1:end);
%
%     % plot accordingly
%     figure(2),
%     patch('Faces',R(ixlo,:),'Vertices',V, ...
%           'FaceColor',graphics_color('Yellow'),'EdgeColor','black');
%     hold on,
%     patch('Faces',R(ixhi,:),'Vertices',V, ...
%           'FaceColor',graphics_color('Sky blue'),'EdgeColor','black');
%
%   See also HES1_SYSTEM, HES1_JACOBIAN.

% S. Engblom 2024-09-30

syms D N M P n
syms alphaD alphaN alphaM alphaP alphan
syms muD muN muM muP mun
syms KM Kn
k = sym('k');
h = sym('h');

f1 = alphaD*n-muD*D;
f2original = alphaN*D-muN*N;
f2 = -muN*N; % (f2original = f2+f2coupling is split into two terms)
f2coupling = alphaN*D;
f3 = alphaM*N/(1+(P/KM)^k)-muM*M;
f4 = alphaP*M-muP*P;
f5 = alphan/(1+(P/Kn)^h)-mun*n;

% RHS function
states = [D N M P n]; % ordering
ode = [f1 f2 f3 f4 f5].';
vars = [D N M P n ...
        alphaD alphaN alphaM alphaP alphan  ...
        muD muN muM muP mun, ...
        KM Kn k h];
fun_ = matlabFunction(ode,'vars',vars);
fun = @(X,alpha,mu,H)fun_(X(1),X(2),X(3),X(4),X(5), ...
                          alpha(1),alpha(2),alpha(3),alpha(4),alpha(5), ...
                          mu(1),mu(2),mu(3),mu(4),mu(5), ...
                          H(1),H(2),H(3),H(4));

% Jacobian function
locJ = jacobian(ode,states);
funJ_ = matlabFunction(locJ,'vars',vars);
funJ = @(X,alpha,mu,H)funJ_(X(1),X(2),X(3),X(4),X(5), ...
                            alpha(1),alpha(2),alpha(3),alpha(4),alpha(5), ...
                            mu(1),mu(2),mu(3),mu(4),mu(5), ...
                            H(1),H(2),H(3),H(4));

% (quite sparse) coupling Jacobian
coupJ = jacobian([0 f2coupling 0 0 0].',states);
funC_ = matlabFunction(coupJ,'vars',vars(numel(states)+1:end));
funC = @(alpha,mu,H)funC_(alpha(1),alpha(2),alpha(3),alpha(4),alpha(5), ...
                          mu(1),mu(2),mu(3),mu(4),mu(5), ...
                          H(1),H(2),H(3),H(4));
