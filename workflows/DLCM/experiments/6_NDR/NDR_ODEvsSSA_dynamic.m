% Notch-delta-reporter model: ODE vs. SSA, growing population.
%
%   The intent with this experiment is to continue on NDR_ODEvsSSA but
%   over a growing population of cells.

% S. Engblom 2018-01-26 (Revision, dynamic case)
% S. Engblom 2018-01-24 (Revision)
% S. Engblom 2017-09-11 (Revision)
% S. Engblom 2016-12-05

if ~exist('Event','var') || ~exist('Etime','var') || ~exist('Etype','var')
  error('Assumes events (Event,Etime,Etype) from a DLCM-run are available.');
end

% sol == 1 or 2 (ODE or SSA)
if ~exist('sol','var'), error('Define variable sol.'); end

% solution recorded at this grid:
tspan = linspace(0,Tend,51);

% delta-notch-reporter DOFs
U = Usave{1};
Notch = spreplace(U,abs(1e-3*param.betaN+1e-4*param.betaN*randn(nnz(U),1)));
Notch = [Notch sparse(size(U,1),1)];
Delta = spreplace(U,abs(1e-3*param.betaD+1e-4*param.betaD*randn(nnz(U),1)));
Delta = [Delta sparse(size(U,1),1)];
Report = sparse(size(U,1),2);
% Format: these are sparse solution vectors over the full mesh
% [P,E,T], where the second column is unused until cells proliferate,
% a filled row then means the corresponding voxel contains two cells
% (which is the maximum).

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Nsave = cell(1,numel(tspan));
Dsave = cell(1,numel(tspan));
Rsave = cell(1,numel(tspan));

% ODE vs. SSA
if sol == 1
  % initial data
  Usave(1) = {U};
  Nsave(1) = {Notch};
  Dsave(1) = {Delta};
  Rsave(1) = {Report};

  % loop over grid in time
  tt = tspan(1);
  k = 1;
  for j = 2:numel(tspan)
    j

    % loop over prerecorded events in tspan(j) <= tt < tspan(j)
    while 1
      % junctional and protrusional contacts
      adof1 = find(U);
      adof2 = find(U == 2);
      nvox = numel(adof1);
      ncell = nvox+numel(adof2);
      % here the convention of contact within the same voxel is very
      % convenient as it ensures that the two off-diagonal blocks are
      % correctly handled
      Nja = [Nj(adof1,adof1) Nj(adof1,adof2); ...
             Nj(adof2,adof1) Nj(adof2,adof2)];
      Npa = [Np(adof1,adof1) Np(adof1,adof2); ...
             Np(adof2,adof1) Np(adof2,adof2)];

      % integrate until next DLCM-layer event
      X0 = full([Notch(adof1,1) Delta(adof1,1) Report(adof1,1); ...
                 Notch(adof2,2) Delta(adof2,2) Report(adof2,2)]);
      odeopts = odeset('Vectorized','on', ...
                       'AbsTol',1e-4,'RelTol',1e-3);
      if k > numel(Etime) || Etime(k) > tspan(j)
        ttend = tspan(j);
      else
        ttend = Etime(k);
      end

      [~,X] = ode23(@NDRrhs,[tt ttend],X0(:),odeopts,Nja,Npa,param);
% $$$       % Jacobian pattern
% $$$       Jp = kron([0 1 0; 1 0 0; 0 1 0],Nja+Npa+speye(ncell)+ ...
% $$$            kron([0 0 0; 0 0 1; 1 0 0],speye(ncell))+ ...
% $$$            speye(3*ncell);
% $$$       odeopts = odeset('Vectorized','on', ...
% $$$                        'AbsTol',1e-4,'RelTol',1e-3,'JPattern',Jp);
% $$$       [~,X] = ode23s(@NDRrhs,[tt Etime(k)],X0(:),odeopts,Nja,Npa,param);

      % update state
      X = reshape(X(end,:),[],3);
      Notch(adof1,1) = X(1:nvox,1);
      Delta(adof1,1) = X(1:nvox,2);
      Report(adof1,1) = X(1:nvox,3);
      Notch(adof2,2) = X(nvox+1:end,1);
      Delta(adof2,2) = X(nvox+1:end,2);
      Report(adof2,2) = X(nvox+1:end,3);
      tt = ttend;

      % exit condition
      if k > numel(Etime) || Etime(k) > tspan(j)
        break;
      end;

      % execute event: update state first
      if Etype(k) == 1
         % movement of a boundary (singly occupied) voxel
        ix = find(Event(:,k) == -1);
        n = find(Event(:,k) == 1);
        % clear empty voxel
        Delta(n,1) = Delta(ix,1); Delta(ix,1) = 0;
        Notch(n,1) = Notch(ix,1); Notch(ix,1) = 0;
        Report(n,1) = Report(ix,1); Report(ix,1) = 0;
      elseif Etype(k) == 2
        % movement of a cell in a doubly occupied voxel
        ix = find(Event(:,k) == -1);
        n = find(Event(:,k) == 1);
        % "movement noise"
        m = 1+(rand > 0.5); % 1 or 2
        if U(n) == 0
          % one cell moved into an empty voxel
          Delta(n,1) = Delta(ix,m);
          Notch(n,1) = Notch(ix,m);
          Report(n,1) = Report(ix,m);
        else
          % one cell moved into an already occupied voxel
          Delta(n,2) = Delta(ix,m);
          Notch(n,2) = Notch(ix,m);
          Report(n,2) = Report(ix,m);
        end
        if m == 1
          Delta(ix,1) = Delta(ix,2);
          Notch(ix,1) = Notch(ix,2);
          Report(ix,1) = Report(ix,2);
        end
        Delta(ix,2) = 0;
        Notch(ix,2) = 0;
        Report(ix,2) = 0;
      elseif Etype(k) == 3
        % proliferation event
        ix = find(Event(:,k) == 1);
         % share Delta-Notch-Reporter between the two cells
         r = rand;
         Delta(ix,2) = r*Delta(ix,1); Delta(ix,1) = (1-r)*Delta(ix,1);
         r = rand;
         Notch(ix,2) = r*Notch(ix,1); Notch(ix,1) = (1-r)*Notch(ix,1);
         r = rand;
         Report(ix,2) = r*Report(ix,1); Report(ix,1) = (1-r)*Report(ix,1);
      end
      % finally, here is the actual event:
      U = U+Event(:,k);
      k = k+1;
    end

    % record values
    Usave(j) = {U};
    Nsave(j) = {Notch};
    Dsave(j) = {Delta};
    Rsave(j) = {Report};
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
  clear umod

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
  umod.seed = 1498;

  % initial data
  Usave(1) = {U};
  Nsave(1) = {Notch};
  Dsave(1) = {Delta};
  Rsave(1) = {Report};

  % loop over grid in time
  tt = tspan(1);
  k = 1;
  for j = 2:numel(tspan)
    j

    % loop over prerecorded events in tspan(j) <= tt < tspan(j)
    while 1
      % junctional and protrusional contacts
      adof1 = find(U);
      adof2 = find(U == 2);
      nvox = numel(adof1);
      ncell = nvox+numel(adof2);
      Nja = [Nj(adof1,adof1) Nj(adof1,adof2); ...
             Nj(adof2,adof1) Nj(adof2,adof2)];
      Npa = [Np(adof1,adof1) Np(adof1,adof2); ...
             Np(adof2,adof1) Np(adof2,adof2)];
      
      % no diffusion
      umod.D = sparse(7*ncell,7*ncell);
      umod.vol = ones(ncell,1);
      umod.sd = ones(ncell,1);
      % *** implementation possibility: use ldata[0,1,2,3] for [N_a N_b D_a
      % D_b] instead

      % determine a suitable micro time-step
      if k > numel(Etime) || Etime(k) > tspan(j)
        ttend = tspan(j);
      else
        ttend = Etime(k);
      end
      dT = ttend-tt;
      X = full([Notch(adof1,1) Delta(adof1,1) Report(adof1,1); ...
                Notch(adof2,2) Delta(adof2,2) Report(adof2,2)]);
      dX = NDRrhs(0,X(:),Nja,Npa,param);
      % allow at most 5% Euler step increase in a norm-wise sense:
      dtau = 0.05*norm(X)./norm(dX);
      % time-step is min of those, but using an integer number of steps:
      dt = min(dtau,dT);
      K = min(ceil(dT/dt),10); % *** to think of
      dt = dT/K;
      umod.tspan = [0 dt];

      % initial conditions
      umod.u0 = [round(VolND*full([Notch(adof1,1); Notch(adof2,2)])'); ...
                 round(VolND*full([Delta(adof1,1); Delta(adof2,2)])'); ...
                 round(VolR*full([Report(adof1,1); Report(adof2,2)])'); ...
                 zeros(4,ncell)];
 
      % take K micro-steps
      for kk = 1:K
        % update junctional and protrusional singalling contacts
        umod.u0(4,:) = umod.u0(1,:)*Nja'-umod.u0(1,:);
        umod.u0(5,:) = umod.u0(1,:)*Npa'-umod.u0(1,:);
        umod.u0(6,:) = umod.u0(2,:)*Nja'-umod.u0(2,:);
        umod.u0(7,:) = umod.u0(2,:)*Npa'-umod.u0(2,:);

        % parse/compile/run SSA-solver:
        umod.seed = umod.seed+1;
        umod = urdme(umod);
        % will not parse/compile again:
        umod.parse = 0;
        umod.compile = 0;
 
        % prepare for next step
        u = reshape(umod.U(:,end),[],ncell);
        umod.u0 = u;
        % (note: umod.u0(4:7,:) is next overwritten)
      end

      % update state
      Notch(adof1,1) = sparse(u(1,1:nvox)/VolND)';
      Delta(adof1,1) = sparse(u(2,1:nvox)/VolND)';
      Report(adof1,1) = sparse(u(3,1:nvox)/VolR)';
      Notch(adof2,2) = sparse(u(1,nvox+1:end)/VolND)';
      Delta(adof2,2) = sparse(u(2,nvox+1:end)/VolND)';
      Report(adof2,2) = sparse(u(3,nvox+1:end)/VolR)';
      tt = ttend;

      % exit condition
      if k > numel(Etime) || Etime(k) > tspan(j)
        break;
      end;

      % execute event: update state first
      if Etype(k) == 1
         % movement of a boundary (singly occupied) voxel
        ix = find(Event(:,k) == -1);
        n = find(Event(:,k) == 1);
        % clear empty voxel
        Delta(n,1) = Delta(ix,1); Delta(ix,1) = 0;
        Notch(n,1) = Notch(ix,1); Notch(ix,1) = 0;
        Report(n,1) = Report(ix,1); Report(ix,1) = 0;
      elseif Etype(k) == 2
        % movement of a cell in a doubly occupied voxel
        ix = find(Event(:,k) == -1);
        n = find(Event(:,k) == 1);
        % "movement noise"
        m = 1+(rand > 0.5); % 1 or 2
        if U(n) == 0
          % one cell moved into an empty voxel
          Delta(n,1) = Delta(ix,m);
          Notch(n,1) = Notch(ix,m);
          Report(n,1) = Report(ix,m);
        else
          % one cell moved into an already occupied voxel
          Delta(n,2) = Delta(ix,m);
          Notch(n,2) = Notch(ix,m);
          Report(n,2) = Report(ix,m);
        end
        if m == 1
          Delta(ix,1) = Delta(ix,2);
          Notch(ix,1) = Notch(ix,2);
          Report(ix,1) = Report(ix,2);
        end
        Delta(ix,2) = 0;
        Notch(ix,2) = 0;
        Report(ix,2) = 0;
      elseif Etype(k) == 3
        % proliferation event
        ix = find(Event(:,k) == 1);
         % share Delta-Notch-Reporter between the two cells
         r = rand;
         Delta(ix,2) = r*Delta(ix,1); Delta(ix,1) = (1-r)*Delta(ix,1);
         r = rand;
         Notch(ix,2) = r*Notch(ix,1); Notch(ix,1) = (1-r)*Notch(ix,1);
         r = rand;
         Report(ix,2) = r*Report(ix,1); Report(ix,1) = (1-r)*Report(ix,1);
      end
      % finally, here is the actual event:
      U = U+Event(:,k);
      k = k+1;
    end

    % record values
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
          sprintf('visual/NDR_ODEvsSSA%d_dynamic.gif',sol), ...
          'delaytime',0.1,'loopcount',0);

return;
