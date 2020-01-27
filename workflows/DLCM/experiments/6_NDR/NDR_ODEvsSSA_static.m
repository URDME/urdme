% Notch-delta-reporter model: ODE vs. SSA.
%
%   The intent with this experiment is to match the ODE-model with an
%   associated stochastic model. The degree of freedom in this
%   interpretation is the system size volume.

% S. Engblom 2018-01-24 (Revision)
% S. Engblom 2017-09-11 (Revision)
% S. Engblom 2016-12-05

% sol == 1 or 2 (ODE or SSA)
if ~exist('sol','var'), error('Define variable sol.'); end

% simulation interval
if ~exist('Tend','var')
  warning('Non-existing variable (Tend). Default used.');
  Tend = 40;
end

% solution recorded at this grid:
tspan = linspace(0,Tend,51);

Nvoxels = size(U,1);
if ~exist('sol_restart','var') || ~sol_restart
  % delta-notch-reporter DOFs
  Notch = spreplace(U,abs(1e-3*param.betaN+1e-4*param.betaN*randn(nnz(U),1)));
  Delta = spreplace(U,abs(1e-3*param.betaD+1e-4*param.betaD*randn(nnz(U),1)));
  Report = sparse(size(U,1),size(U,2));
else
  % alternative, start from previous:
  clear sol_restart % force explicit reuse
  Notch = Nsave{end};
  Delta = Dsave{end};
  Report = Rsave{end};
end

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Nsave = cell(1,numel(tspan));
Dsave = cell(1,numel(tspan));
Rsave = cell(1,numel(tspan));

% ODE vs. SSA
if sol == 1
  % Jacobian pattern
  Jp = kron([0 1 0; 1 0 0; 0 1 0],Nj+Np+speye(Nvoxels))+ ...
       kron([0 0 0; 0 0 1; 1 0 0],speye(Nvoxels))+ ...
       speye(3*Nvoxels);
  odeopts = odeset('OutputFcn',@report,'Vectorized','on', ...
                   'AbsTol',1e-4,'RelTol',1e-3,'JPattern',Jp);  
  X0 = full([Notch Delta Report]);
  [~,X] = ode15s(@NDRrhs,tspan,X0(:),odeopts,Nj,Np,param);

  % record values
  X = reshape(X',[],3,size(X,1));
  for j = 1:numel(tspan)
    Usave(j) = {U};
    Nsave(j) = {sparse(X(:,1,j))};
    Dsave(j) = {sparse(X(:,2,j))};
    Rsave(j) = {sparse(X(:,3,j))};
  end
else
  % setup URDME model in variables [N D R N_a N_b D_a D_b]

  % the scaling variables for (N,D) and R enables a comparision to the
  % deterministic model
  if ~exist('VolND','var') || ~exist('VolR','var')
    error('Define volume variables.');
  end
  %   N' = betaN-<D_in>*N/kt-D*N/kc-N
  %   D' = betaD*G(R)-D*<N_in>/kt-D*N/kc-D
  %   R' = betaR*F(<D_out>*N)-R
  %
  % scaling through by [VolND VolR] yields,
  %
  %   (VolND*N)' = VolND*betaN-<D_in>*VolND*N/kt-D*VolND*N/kc-VolND*N
  %   (VolND*D)' = VolND*betaD*G(R)-VolND*D*<N_in>/kt-
  %                  VolND*D*N/kc-VolND*D
  %   (VolR*R)' = VolR*betaR*F(<D_out>*N)-VolR*R
  %
  % where
  %
  %   G(x) = 1/(1+x^m), (m = 2),
  %   F(x) = x^s/(kRS+x^s), (s = 2).
  %
  % And so, in terms of the scaled variables [N] = VolND*N, [D] =
  % VolND*D, [R] = VolR*R,
  %
  %   [N]' = VolND*betaN-<[D]_in>*[N]/(kt*VolND)-[D]*[N]/(kc*VolND)-[N]
  %   [D]' = VolND*betaD*G([R]/VolR)-[D]*<[N]_in>/(kt*VolND)-
  %            [D]*[N]/(kc*VolND)-[D]
  %   [R]' = VolR*betaR*F(<[D]_out>*[N]/VolND^2)-[R]
  
  % generate URDME model

  % no diffusion
  clear umod
  umod.D = sparse(7*Nvoxels,7*Nvoxels);
  umod.vol = ones(Nvoxels,1);
  umod.sd = ones(Nvoxels,1);
  % *** implementation possibility: use ldata[0,1,2,3] for [N_a N_b D_a
  % D_b] instead

  % kinetic model
  r = cell(1,9);
  % (1) birth N, (2) degradation N (trans), (3) mutual degradation
  % (N,D), (4) degradation N
  r{1} = '@ > betaN*vol > N';
  r{2} = 'N > N*(wa*D_a+wb*D_b)/kt > @';
  r{3} = 'N+D > N*D/(kc*vol) > @';
  r{4} = 'N > N > @';

  % (5) birth D, downregulated by R, (6) degradation D (trans), (7)
  % degradation D
  r{5} = '@ > betaD*vol/(1.0+pow(R/(VolR*vol),2.0)) > D';
  r{6} = 'D > D*(wa*N_a+wb*N_b)/kt > @';
  r{7} = 'D > D > @';

  % (8) birth R, upregulated by (trans-)D and N, (9) death R
  r{8} = '@ > betaR*vol/(kRS/pow((qa*D_a+qb*D_b)/VolND*(N/(vol*VolND)),2.0)+1.0) > R';
  r{9} = 'R > R > @';

  spec = {'N' 'D' 'R' 'N_a' 'N_b' 'D_a' 'D_b'};
  rate = {'VolND' VolND 'VolR' VolR ...
          'betaN' VolND*param.betaN 'betaD' VolND*param.betaD ...
          'betaR' VolR*param.betaR ...
          'kt' VolND*param.kt 'kc' VolND*param.kc 'kRS' param.kRS ...
          'wa' param.wa 'wb' param.wb 'qa' param.qa 'qb' param.qb};
  [~,umod.N,umod.G] = rparse(r,spec,rate,'source/NDR_ODEvsSSA.c');
  umod.solver = 'ssa';
  umod.propensities = 'source/NDR_ODEvsSSA.c';
  umod.seed = 1709;

  % initial data
  Usave(1) = {U};
  Nsave(1) = {Notch};
  Dsave(1) = {Delta};
  Rsave(1) = {Report};

  for j = 2:numel(tspan)
    j

    % determine a suitable micro time-step
    dT = tspan(j)-tspan(j-1);
    X = full([Notch; Delta; Report]);
    dX = NDRrhs(0,X,Nj,Np,param);
    % allow at most 5% Euler step increase in a norm-wise sense:
    dtau = 0.05*norm(X)./norm(dX);
    % time-step is min of those, but using an integer number of steps:
    dt = min(dtau,dT);
    K = min(ceil(dT/dt),10); % *** to think of
    dt = dT/K;
    umod.tspan = [0 dt];

    % initial conditions
    umod.u0 = [round(VolND*full(Notch)'); ...
               round(VolND*full(Delta)'); ...
               round(VolR*full(Report)'); ...
               zeros(4,Nvoxels)];
 
    % take K micro-steps
    for k = 1:K
      % update junctional and protrusional singalling contacts
      umod.u0(4,:) = umod.u0(1,:)*Nj'-umod.u0(1,:);
      umod.u0(5,:) = umod.u0(1,:)*Np'-umod.u0(1,:);
      umod.u0(6,:) = umod.u0(2,:)*Nj'-umod.u0(2,:);
      umod.u0(7,:) = umod.u0(2,:)*Np'-umod.u0(2,:);

      % parse/compile/run SSA-solver:
      umod.seed = umod.seed+1;
      umod = urdme(umod);
      % will not parse/compile again:
      umod.parse = 0;
      umod.compile = 0;
 
      % prepare for next step
      u = reshape(umod.U(:,end),[],Nvoxels);
      umod.u0 = u;
      % (note: umod.u0(4:7,:) is next overwritten)
    end

    % solution at grid point is now available
    Notch = sparse(u(1,:)/VolND)';
    Delta = sparse(u(2,:)/VolND)';
    Report = sparse(u(3,:)/VolR)';
    Usave(j) = {U};
    Nsave(j) = {Notch};
    Dsave(j) = {Delta};
    Rsave(j) = {Report};
  end
end

% display final frame
figure(sol), clf,
i = numel(Usave);
NDRplot(Dsave{i},Usave{i},V,R);

return;

% movie visualization
mx = full(max(Dsave{end}(:)));
figure(sol),
M = struct('cdata',{},'colormap',{});
for i = 1:numel(Usave)
  clf,
  NDRplot(Dsave{i},Usave{i},V,R,mx);
  M(i) = getframe(gcf);
end

movie2gif(M,{M([1:2 end]).cdata}, ...
          sprintf('visual/NDR_ODEvsSSA%d.gif',sol), ...
          'delaytime',0.1,'loopcount',0);

return;
