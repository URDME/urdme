function U = mexdlcm(mexhash,tspan,u0,D,N,G,vol,ldata,gdata, ...
                    data_time,ldata_time,gdata_time, ...
                    sd,reportl,seed,K,I,S,solverargs)
% (for help, type 'help dlcm')

% E. Blom 2024-04-25 (continuation from DLCM workflow)

%% (0) Unpack input
if ~isscalar(seed)
  error('Solver cannot handle replicas')
end
rng(seed);

% load necessary parameters from ldata and gdata
P = ldata(1:2,:);           % map: voxel to position
nquants = gdata(1);         % number of micro-quantities
ntypes = gdata(2);          % number of cell types

opts = struct(solverargs{:});     % cell to struct
gradquotient = opts.gradquotient; % voxel edge to neighbour distance ratio

% create function handles for evaluating reaction rates and Q rhs:s
mexrhs = str2func([opts.mexname '_mexrhs']);

Rates = opts.Rates;                 % migration event rates
Drate = opts.Drate;                 % migration rate scaling
D = opts.D;                         % Laplacian
Ne = opts.Ne;                       % Neighbouring matrix

%% (0a) Internal states
% number of internal states (including cell voxel and type)
ninternal = 2;                      % states added if mumod exists...

% load micro-state dynamics, if possible (as URDME struct)
internal_rates_on = isfield(opts.mumod, 'N'); % existence of N -> inner SSA
internal_is_discrete = true;        % discrete internal states (default)
if ~isempty(opts.mumod)
  mumod = opts.mumod;
  ninternal = ninternal + (size(mumod.u0,1))/2;
  if strcmp(opts.internal_state,'cont')   % continuous internal states
    % compile once to be able to use mexrhs for internal events with mumod
    mumod_tmp = urdme(mumod,'solver','uds', ...
         'solve',0);
    mexrhs_internal = str2func([mumod_tmp.mexname '_mexrhs']);
    internal_is_discrete = false;   % continuous internal states
  end
end

% load additional internal functions (empty ones are not used)
maxdt_fun = opts.maxdt_fun; % function for SSA max. time-interval
ldata_fun = opts.ldata_fun; % function for ldata in SSA solver

%% (0b) Curvature & Surface tension
% load elliptic projection operators, if possible
if ~isempty(opts.curv)
  curv = struct(opts.curv{:});
  curvop_exists = true;
  SLa.X = curv.SLa;
  M     = curv.M;
  Ax    = curv.Ax;
  Ay    = curv.Ay;
  is3D  = false;
  if isfield('Az',curv) % if 3D
    Az  = curv.Az;
    is3D = true;
  end
  sigma = curv.sigma;
else
  curvop_exists = false;        % skip curvature and Y-L evaluation
end

%% (0c) Check inputs
% check N -- should correspond to single cell switching type or
% proliferating or death
[~, j] = find(N);  % find nnz occurences per column
if ~all(abs(N(:))<=2)...              % no element larger than size two
        || ~all(N(:)>=-1)...          % no element smaller than minus one
        || ~all(groupcounts(j) <= 2)  % no column has more than 2 nz rows
  error('Stoichiometric matrix must correspond to single-cell event.')
end

%% (0d) Factorize operators
% explicitly get extdof (boundary nodes) from sd
extdof = find(sd==0);

% factorize Laplacian
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
La.X = D; % use the negative Laplacian
Lai = fsparse(extdof,extdof,1,size(La.X));
La.X = La.X-Lai*La.X+Lai;
[La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

% LU-factorize elliptic projection operations
if curvop_exists
  [SLa.L,SLa.U,SLa.p,SLa.q,SLa.R] = lu(SLa.X,'vector');
end

%% (0e) Allocate save states & define cell population and internal states
% representation of cell population, internal states, etc.
nsavestates = ntypes + 2 + (ninternal-2)*2; % cell types, states, gl. index
Usave = zeros(nsavestates,numel(u0(1,:)),numel(tspan));
% ... and similarly for global quantities and cell index-voxel map
Qsave = zeros(nquants,numel(u0(1,:)),numel(tspan));

Tend = tspan(end);          % Simulation interval
tt = tspan(1);

% U is an Nvoxels-by-1 array representing the cell population during
% simulation
U = sum(u0(1:ntypes,:),1);
if any(U>2)
  error('Initial population voxel capacity breached ')
end
i = 1;                      % index used to save states at proper times
neigh = full(sum(Ne,2));    % nr of neighbors to each voxel

% Cell indexing and data. IX is an Nvoxels-by-Ncapacity array containing
% each cell indices: cell # in voxel i has unique index IX(i, #).
Nvoxels = numel(U);
IX = ones(Nvoxels, 2);              % empty voxels have index=1
adof = find(U>0);
sdof = find(U>1);
IX(adof, 1) = 2:(numel(adof)+1);            % give index to each #=1
IX(sdof, 2) = numel(adof)+1+(1:numel(sdof));% give index to each #=2
ix_max = max(IX(:));                        % keeps track of indices used

% Cell data X is a 2+ninternal-by-Nvoxels*Ncapacity array. X(n, idx)
% contains the n:th internal data of cell with index idx. I.e, cell # in
% voxel that has idx = IX(i, #). 1st and 2nd internal data is reserved for
% which voxel the cell is in and which cell type it is, respectively.
X = zeros(ninternal, Nvoxels*2);   % allocate for max nr of cells
X(1:2, 1) = [0 0];                         % empty voxel 'data'
X(1, IX(adof, 1)) = adof;                  % which voxel cell is in
X(1, IX(sdof, 2)) = sdof;
for n = 1:ntypes                           % assign cell types
  X(2, IX(u0(n,:)>0, 1)) = n; % singly occ.
  X(2, IX(u0(n,:)>1, 2)) = n; % double occ.
end

% Load internal cell states into X for n = 3, ..., ninternal
for n = 3:ninternal
  X(n, IX(:,1)) = mumod.u0(2*(n-2)-1, :);  % 1st cell
  X(n, IX(:,2)) = mumod.u0(2*(n-2), :);    % 2nd cell
  Usave(ntypes+2*(n-2)-1,:,1) = X(n, IX(:,1));
  Usave(ntypes+2*(n-2),:,1) = X(n, IX(:,2));
end

% Q is an Nvoxels-by-nquant array containing each continuous quantity
% modeled by the Qrhs function in respective column. Q(:,1) is reserved
% for pressure
Q = u0(ntypes+(1:nquants), :)';

% save initial cell data...
for n = 1:ntypes       % save voxel data: number of type n in each voxel
  Usave(n,:,1) = u0(n,:);
end
Usave(nsavestates-1, :, 1) = IX(:,1);       % global cell index, 1st cell
Usave(nsavestates, :, 1) = IX(:,2);         % 2nd cell
% ... and save intitial global quantities
for q = 1:nquants                           % save global quantities
  Qsave(q,:,1) = Q(:,q);
end

% ...after this u0 is no longer needed.
clear u0;

%% (1) solve
U = U';    % urdme parser and dlcm solver use different formats for U...
while tt < tspan(end)

  %% (1a) classify the Degrees-Of-Freedoms (DOFs for short)
  % Find all voxels with living cells and classify voxels from which
  % movement is allowed, bdof_m and sdof_m.

  % active DOFs: occupied voxels
  adof = find(U);

  % boundary DOFs _which may move_: containing one cell per voxel and
  % with an empty voxel besides it
  bdof_m = find(Ne*(U > 0) < neigh & U == 1);
  % also consider those singly occupied voxels neighbouring other cell
  % types as bdof_ms
  for n = 1:ntypes
    bdof_mn = find(Ne*(X(2, IX(:,1)) == n)' < neigh & ...
        (X(2, IX(:,1)) == n)' & U == 1);
    bdof_m = union(bdof_m, bdof_mn);
  end
  % required for the special case when all size(bdof_mn) = 1:
  bdof_m = reshape(bdof_m, numel(bdof_m), 1);

  % source DOFs _which may move_: with a voxel containing less number of
  % cells next to it
  sdof_m = find(Ne*(U > 1) < neigh & U > 1);

  % injection DOFs: the layer of voxels "just outside" adof, used for
  % population BCs.
  idof = find(Ne*(U ~= 0) > 0 & U == 0);

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];
  % (after this adof is not used)

  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';
  [bdof_m_,sdof_m_,idof_] = ...
      map(Adof_,Adof,bdof_m,sdof_m,idof);

  % get QI -- an Ncells-by-Ncapacity-by-ninternal array with every
  % cell's internal state.
  QI = zeros(numel(U),2,ninternal);
  QI(:,1,:) = X(1:ninternal, IX(:, 1))';
  QI(:,2,:) = X(1:ninternal, IX(:, 2))';

  %% (1b) Evaluate curvature per cell type

  % Evaluate curvature by three steps of elliptic projection. Assuming that
  % voxels are homotypic, we have here that tdof_(type==n) gives dofs of
  % type n in the local enumeration.
  [type, tdof_] = find(X(2, IX(Adof,1)) == (1:ntypes)');
  curvature = zeros(numel(U), 1);
  if curvop_exists          % (else curvature evaluation is disabled)
    for ti = 1:ntypes % Get curvatures separately for each cell type
      Utype = zeros(numel(U), 1);
      Utype(Adof(tdof_(type==ti))) = 1;
      % 1/3: Smooth cell population...
      S(SLa.q,1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\(M*Utype)));
      % 2/3: ...get its gradient...
      Sx(SLa.q,1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Ax*S));
      Sy(SLa.q,1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Ay*S));
      if is3D
        Sz(SLa.q,1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Az*S));
        nrm = 1./sqrt(Sx.^2 + Sy.^2 + Sz.^2);
        Sz_norm = -Sz.*nrm;
      else
        nrm = 1./sqrt(Sx.^2 + Sy.^2);  % normalize
      end
      Sx_norm = -Sx.*nrm;
      Sy_norm = -Sy.*nrm;
      % 3/3: ...get its divergence = curvature
      curv1(SLa.q, 1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Ay*Sy_norm));
      curv2(SLa.q, 1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Ax*Sx_norm));
      if is3D
        curv3(SLa.q, 1) = SLa.U\(SLa.L\(SLa.R(:,SLa.p)\Az*Sz_norm));
        % keep only curvature over voxels with cell type ti...
        curvature(Adof(tdof_(type==ti)),1) = curv3(Adof(tdof_(type==ti))) + ...
          curv1(Adof(tdof_(type==ti))) + curv2(Adof(tdof_(type==ti)));
      else
        % keep only curvature over voxels with cell type ti...
        curvature(Adof(tdof_(type==ti)),1) =  ...
          curv1(Adof(tdof_(type==ti))) + curv2(Adof(tdof_(type==ti)));
      end
      % ...and its idofs (prioritising larger type values at overlaps)
      [a, ~] = find(Ne(idof,Adof(tdof_(type==ti))));
      curvature(idof(a),1) = curv1(idof(a)) + curv2(idof(a));
    end
  end

  %% (1c) Calculate micro-environment quantities
  % Starting with pressure and the Y-L conditions

  % get all rates as follows (including quantity rhs:s and bc:s):
  y = (X(2, IX(:, 1))'==1:ntypes)';   % current state
  y = y + (X(2, IX(:, 2))'==1:ntypes)';
  y = [y; Q'];                        % previous micro-quantities
  % Use ldata for QI
  for l = 3:ninternal     % (ldata 1,2 holds position map)
    ldata(l,:) = QI(:,1,l);               % 1nd cell
    ldata(l+ninternal-2,:) = QI(:,2,l);   % 2nd cell
  end
  R = mexrhs(mexhash,tt,y,size(N,2), ...
           vol,ldata,gdata, ...
           ldata_time,gdata_time,sd, ...
           [],[],[]);
  R = reshape(R,[],numel(vol)); % (Mreactions-by-Nvoxels)

  % all Q:s have exactly one R associated to them:
  source = R(end-nquants+1,Adof);    % pressure sources (first defined Q)

  % Construct Young-Laplace (Y-L) pressure drop conditions between all
  % voxels with cells of different type (local enumeration). The Y-L
  % conditions are formulated as Lagrange conditions to the pressure
  % system of eqs. Letting ti-->tj indicate cell type ti connected to cell
  % type tj, it necessary that ti>tj to get the appropriate curvature sign.
  Lag = [];                             % Lagrange conditions for Y-L
  Lag_rhs = [];                         % Lagrangre system rhs
  used_idx = [];                        % - to ensure unique nodes
  if curvop_exists                      % (skip if surface tension absent)
  for ti = ntypes:-1:1                  % from type ti...
    tdofi_ = tdof_(type==ti);
    for tj = (ti-1):-1:1                % ... to type tj, with ti>tj
      tdofj_ = tdof_(type==tj);
      % find interface dofs between cells of different type
      [i_, j_] = find(Ne(Adof(tdofi_), Adof(tdofj_)));
      in_ = tdofi_(i_);
      % ensure unique in_:s, to avoid singularity ...
      [in_, keep] = setdiff(in_, used_idx, 'stable');
      jn_ = tdofj_(j_(keep));
      used_idx = union(used_idx,in_); % ...and remember used in_:s
      if ~isempty(in_)                % if valid interface exists:
        % use mean curvature between dofs (removes discretization issue)
        C = -0.5*curvature(Adof(in_)) + 0.5*curvature(Adof(jn_));
        % append Lagrange conditions
        Lag = [Lag; fsparse((1:numel(in_))', ...
          [in_ jn_],[-1 1],[numel(in_) numel(Adof)])];
        Lag_rhs = [Lag_rhs; sigma(ti+1,tj+1)*C];
      end
    end
    % consider single voxels for the boundary of the full population
    tdofi_ = tdof_(type==ti);
    % find idofs to cell type ti...
    [in, ~] = find(Ne(idof, Adof(tdofi_)));
    % ... and set BC there:
    source(idof_(in)) = sigma(ti+1,1)*curvature(idof(in));
    % Note: smaller type value overrides on duplicate idofs!
  end
  end
  % *** DEBUG: use to plot dof connections
  %$$$ plot([P(1,tdof1(i1)); P(1,tdof2(j1))], [P(2,tdof1(i1)); ...
  %$$$ P(2,tdof2(j1))], 'ko-')

  % Impose outer BC (against 'emptiness') on tumor boundary directly
  L = D(Adof,Adof);
  Lai = fsparse(idof_,idof_,1,size(L));
  L = L-Lai*L+Lai;

  % Add Lagrange multipliers to system of eqs
  A = [[L Lag']; ...
    [Lag sparse(size(Lag,1),size(Lag,1))]];
  b = [source(:); Lag_rhs];
  Pr = A\b;                           % solve...
  Q(Adof,1) = Pr(1:numel(Adof),1);    % ...keep only relevant values

  % Calculate other Laplace quantities (RHS:s evaluated at current state)
  for q = 2:nquants
    source = R(end-nquants+q,:); % get rhs (bcs are values at sd == 0)
    Q(La.q, q) = La.U\(La.L\(La.R(:,La.p)\source')); % solve
  end

  %% (1d) Calculate event rates
  % Find the movement rates from bdof_m and sdof_m, and reaction event
  % rates defined by user input in Rates and stoichiometric  matrix N.

  % get reaction rates (again; necessary to ensure proper dependence on Q)
  % (ldata not yet changed before last update)
  y(end-nquants+1:end,:) = Q';        % get current micro-quantities
  R = mexrhs(mexhash,tt,y,size(N,2), ...
           vol,ldata,gdata, ...
           ldata_time,gdata_time,sd, ...
           [],[],[]);
  R = reshape(R,[],numel(vol)); % (Mreactions-by-Nvoxels)

  % Get all migration event rates (local enumeration)
  cond = Rates(U(Adof),Q(Adof,:),QI(Adof,:,:),P(:,Adof),tt);

  % Cell movement from voxel with U=1
  % save iib, jjb_ to find which cell moves where in event execution
  [iib, jjb_] = find(Ne(bdof_m,Adof));  % neighbors
  % keep only movement between cells of same type or into void,
  % homotypic voxels assumed
  keep = find(X(2, IX(bdof_m(iib),1)) == X(2, IX(Adof(jjb_),1)) ...
            | X(2, IX(Adof(jjb_),1)) == 0);
  iib = reshape(iib(keep),[],1); jjb_ = reshape(jjb_(keep),[],1);

  % Cell movement from voxel with U=2
  % save iis, jjs_ to find which cell moves where in event execution
  [iis, jjs_] = find(Ne(sdof_m,Adof));  % neighbors
  % keep only movement between cells of same type or into void,
  % homotypic voxels assumed
  keep = find(X(2, IX(sdof_m(iis),1)) == X(2, IX(Adof(jjs_),1)) ...
            | X(2, IX(Adof(jjs_),1)) == 0);
  iis = reshape(iis(keep),[],1); jjs_ = reshape(jjs_(keep),[],1);

  % Movement rates defined by the user input migration pressure
  moveb = zeros(numel(iib),1);                  % migration rates, singly
  moves = zeros(numel(iis),1);                  % doubly occupied
  Drb = Drate(U(bdof_m(iib)),U(Adof(jjb_)), ... % Drates for bdof_m
            Q(bdof_m(iib),:),QI(bdof_m(iib),:,:),P(:,bdof_m(iib)),tt);
  Drs = Drate(U(sdof_m(iis)),U(Adof(jjs_)), ... % Drates for sdof_m
            Q(sdof_m(iis),:),QI(sdof_m(iis),:,:),P(:,sdof_m(iis)),tt);
  for r = 1:numel(cond)         % loop over all migration rates
    Pr = cond{r};               % migration pressure
    if isscalar(gradquotient)   % structured mesh
      moveb = moveb + max(Drb{r}*gradquotient...
          .*(Pr(bdof_m_(iib)) - Pr(jjb_)), 0);
      moves = moves + max(Drs{r}*gradquotient...
          .*(Pr(sdof_m_(iis)) - Pr(jjs_)), 0);
    else                        % unstructured mesh
      moveb = moveb + max(Drb{r}.*diag(gradquotient(bdof_m(iib), ...
          Adof(jjb_))).*(Pr(bdof_m_(iib)) - Pr(jjb_)), 0);
      moves = moves + max(Drs{r}.*diag(gradquotient(sdof_m(iis), ...
          Adof(jjs_))).*(Pr(sdof_m_(iis)) - Pr(jjs_)), 0);
    end
  end
  % *** DEBUG: use to plot resulting pressure gradient arrows:
  %$$$ quiver(P(1,bdof_m(iib)), P(2,bdof_m(iib)), moveb'.*(P(1,Adof(jjb_)) ...
  %$$$   - P(1,bdof_m(iib))), moveb'.*(P(2,Adof(jjb_)) - P(2,bdof_m(iib))))

  % sum of all reaction events/voxel - find which event later
  reaction = sum(R(1:end-nquants,Adof),1)';

  intens = [moveb; moves; reaction];

  % waiting time until next event and the event itself
  lambda = sum(intens);
  dt = -reallog(rand)/lambda;
  rnd = rand*lambda;
  cum = cumsum(intens);
  ix_ = find(rnd < cum, 1, 'first');
  % (now ix_ points to the intensity which fired first;

  %% (1e) Update cell internal states
  % Update the internal state within the Gillespie time interval dt

  Nsteps = 1;     % ensure that data is always saved in proper intervals
  ddt = dt;       % -||-

  % Find micro-time step ddt, dynamic or as fixed user input
  if internal_rates_on
    mumod.sd = zeros(1,numel(U));
    mumod.sd(U>0) = 1;                   % only allow reactions on sd == 1
    mumod.sd(U==2) = 2;
    % Get the complex signals as ldata (assumed constant during ddt)
    ldata_f = ldata_fun(U,Q,QI,Ne);
    for l = 1:numel(ldata_f)
      mumod.ldata(l,:) = ldata_f{l};
    end
    maxdt = maxdt_fun(U,Q,QI,ldata_f);
    % IC for the SSA over this interval
    for n = 3:ninternal
      mumod.u0(2*(n-2)-1,:) = QI(:,1,n);      % 1st cell
      mumod.u0(2*(n-2),:) = QI(:,2,n);        % 2nd cell
    end

    % Sub-time step, ddt
    ddt = min(Tend-tt,dt);
    Nsteps = max(ceil(ddt/maxdt),1);  % even with maxdt > ddt inf, run once
    ddt = ddt/Nsteps;
    mumod.tspan = [0 ddt];            % SSA solver timespan
  end

  % Update internal states (always runs at least once/dt to save states)
  for subt = 1:Nsteps

    % record values as needed
    if tspan(i+1) <= tt+ddt
      iend = i+find(tspan(i+1:end) <= tt+ddt,1,'last');
      % save cell states
      dim = ones(1,numel(i+1:iend));  % ensure correct dimensions in save
      for n = 1:ntypes                % save local cell data
        Usave(n, :, i+1:iend) = ...
          ((X(2, IX(:,1))==n) + (X(2, IX(:,2))==n))'*dim;
      end
      for n = 3:ninternal
        Usave(ntypes+2*(n-2)-1,:,i+1:iend) = X(n, IX(:,1))'*dim;
        Usave(ntypes+2*(n-2),:,i+1:iend) = X(n, IX(:,2))'*dim;
      end
      Usave(nsavestates-1, :, i+1:iend) = IX(:,1)*dim; % global cell indices
      Usave(nsavestates, :, i+1:iend) = IX(:,2)*dim;

      % save global quantities
      for q = 1:nquants                     % save global quantities
        Qsave(q,:,i+1:iend) = Q(:,q)*dim;
      end

      i = iend;
    end

    % check again, since we must allow above state-save in between
    if internal_rates_on
      % Update internal states
      if internal_is_discrete
        % Discrete events; using URDMEs SSA solver
        mumod = urdme(mumod,'solve',1,'parse',0,'compile',0,'solver','ssa');
        mumod.seed = mumod.seed + 1;

        % reshape output to be the new input
        u = reshape(mumod.U(:,end),size(mumod.u0,1),size(mumod.u0,2));
        mumod.u0 = u;                                 % update u0...
      else
        % Continuous events; Euler forward using mexrhs_internal rates
        R_internal = mexrhs_internal(mumod.mexhash,tt,mumod.u0, ...
                size(mumod.N,2),mumod.vol,mumod.ldata,mumod.gdata, ...
               mumod.ldata_time,mumod.gdata_time,mumod.sd, ...
               [],[],[]);
        R_internal = reshape(R_internal,[],numel(mumod.vol));
        mumod.u0 = mumod.u0 + ddt*mumod.N*R_internal; % update u0...
      end

      % ...and copy internal u0 to corresponding 'dlcm'-states
      for n = 3:ninternal
        X(n, IX(:,1)) = mumod.u0(2*(n-2)-1, :);  % 1st cell
        X(n, IX(:,2)) = mumod.u0(2*(n-2), :);    % 2nd cell
      end

      % update QI here for easier input to ldata_fun
      QI(:,1,3:ninternal) = [X(3:ninternal, IX(:, 1))'];
      QI(:,2,3:ninternal) = [X(3:ninternal, IX(:, 2))'];

      % update internal ldata (the more involved user-defined signals)
      ldata_f = ldata_fun(U,Q,QI,Ne);
      for l = 1:numel(ldata_f)
        mumod.ldata(l,:) = ldata_f{l};
      end
    end

    tt = tt + ddt;
  end

  % equilibrium reached at tt, no more events to execute
  if dt == inf
    warning('equilibrium reached: dt = infinity')
    break;
  end

  %% (1f) Execute discrete event
  % Find whether Gillespie event was migration or reaction and execute it.
  % Reactions events split into phenotype switch or birth event.

  % If migration event:
  if ix_ <= numel(moveb) + numel(moves)
    if ix_ <= numel(moveb) % 1-->0
      % find which cell moves and to which voxel it moves
      ix = bdof_m(iib(ix_));              % voxel to move from ...
      n = Adof(jjb_(ix_));                % voxel to move to
    elseif ix_ <= numel(moveb) + numel(moves) % 2-->0 or 1
      % find which cell moves and to which voxel it moves
      ix_ = ix_ - numel(moveb);           % index of ...
      ix = sdof_m(iis(ix_));              % ...voxel to move from
      n = Adof(jjs_(ix_));                % voxel to move to
    end
    % execute MOVE event: move from ix to n
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;

    % move cell index accordingly (moving from 1st position as a FIFO queue)
    IX(n, U(n)) = IX(ix, 1);          % cell moved from ix-->n
    X(1,IX(n, U(n))) = n;             % update which voxel the cell is in
    IX(ix, 1) = IX(ix, 2);            % move cells down "queue" in voxel
    IX(ix, U(ix)+1) = 1;              % no cell left at top position
  % If reaction event:
  else
    % phenotype switch
    ix_ = ix_ - numel(moveb) - numel(moves);
    ix = Adof(ix_);
    Nr = N(1:ntypes,1:end-nquants);  % assumes reactions are defined first!
    % find which cell type to switch to
    rates = R(1:end-nquants, ix);       % all reaction rates at ix
    event = find(cumsum(rates) > rand*sum(rates),1,'first');
    to = find(Nr(:,event) > 0);         % remember, in case of switch event
    event = sum(Nr(:,event));           % -1 :death; 0: switch; 1: birth
    if event == -1      % death event
      X(:, IX(ix, U(ix))) = 0;          % ...clean cell data
      IX(ix, U(ix)) = 1;                % ...deallocate indices
      U(ix) = U(ix) - 1;
    else                % phenotype switch
      % change type of all cells in voxel (preserves homotypicity):
      % and do even if it's a birth event to allow U1 --> U2+U2 events
      % (standard birth events simply switches type to type it already has)
      X(2, IX(ix, 1:U(ix))) = to;       % other data remain same
    end
    if event == 1       % birth event (always preceded by switch event)
      % add cell
      U(ix) = U(ix) + 1;
      % use new unique index
      ix_max = ix_max + 1;
      IX(ix, U(ix)) = ix_max;
      X(:,IX(ix, U(ix))) = X(:,IX(ix, U(ix)-1)); % copy cell data
      tmp_ix = sort(IX(Adof));
      if sum(diff(tmp_ix(tmp_ix>1))==0)
        error('bug: non-unique cell indices!')
      end
    end
  end
  if ~all(U<=2)
    error('Voxel cell capacity limit breached!')
  end
end

%% (2) Save data
% U is ntypes+nquants+2*(ninternal-2)+1-by-ncells-by-ntimestamps. The first
% ntypes rows hold the number of cells per cell type in each voxel, the
% next nquants rows hold micro-environment quantity values per voxel, the
% next 2*(ninternal-2) rows hold the internal states for two cells per
% voxel, excluding cell voxel index and type (they are implied by other
% data). All data is saved per timestamp in tspan.
U = Usave(1:ntypes,:,:);
U(end+1:end+nquants, :,:) = Qsave;
U(end+1:end+(ninternal-2)*2,:,:) = Usave(ntypes+1:ntypes+(ninternal-2)*2,:,:);
U(end+1:end+2,:,:) = Usave(ntypes+1:ntypes+2,:,:) ;% saved in private.dlcm
