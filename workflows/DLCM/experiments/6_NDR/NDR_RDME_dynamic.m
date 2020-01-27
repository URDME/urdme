% Notch-delta-reporter model: RDME model, growing population.

% S. Engblom 2018-01-26 (Revision, dynamic case)
% S. Engblom 2018-01-24 (Revision)
% S. Engblom 2017-09-12 (RDME)

if ~exist('Event','var') || ~exist('Etime','var') || ~exist('Etype','var')
  error('Assumes events (Event,Etime,Etype) from a DLCM-run are available.');
end

% solution recorded at this grid:
tspan = linspace(0,Tend,51);

% delta-notch-reporter DOFs
U = Usave{1};
Notch = spreplace(U,abs(1e-3*param.betaN+1e-4*param.betaN*randn(nnz(U),1)));
Notch = [Notch sparse(size(U,1),1)];
Delta = spreplace(U,abs(1e-3*param.betaD+1e-4*param.betaD*randn(nnz(U),1)));
Delta = [Delta sparse(size(U,1),1)];
Report = sparse(size(U,1),2);
ID = [spreplace(U,1:nnz(U)) sparse(size(U,1),1)];
% Format: these are sparse solution vectors over the full mesh
% [P,E,T], where the second column is unused until cells proliferate,
% a filled row then means the corresponding voxel contains two cells
% (which is the maximum).

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Nsave = cell(1,numel(tspan));
Dsave = cell(1,numel(tspan));
Rsave = cell(1,numel(tspan));

% RDME-layer description
[P_,E_,T_] = RDMElayer;

% assemble the diffusion part

% setup URDME model in variables [N D R N_a N_b D_a D_b]
clear umod
Diff = 1;
umod = pde2urdme(P_,T_,{Diff Diff Diff 0 0 0 0});
% (note: variables [N_a N_b D_a D_b] are signal variables and are in a
% sense 0-dimensional)
% *** implementation possibility: use ldata[0,1,2,3] instead

% inner layer mesh
umod.sd = floor(umod.sd);
Nvoxels_ = numel(umod.sd);

% scale so that each cell is unit total volume
umod.vol = umod.vol/sum(umod.vol);

% the scaling variables for (N,D) and R enables a comparision to the
% deterministic model
if ~exist('VolND','var') || ~exist('VolR','var')
  error('Define volume variables.');
end

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
[~,umod.N,umod.G] = rparse(r,spec,rate,'source/NDR_RDME.c');
umod.solver = 'nsm';
umod.propensities = 'source/NDR_RDME.c';
umod.seed = 1909;
umod.report = 0; % ***

% initial data - outer layer
Usave(1) = {U};
Nsave(1) = {Notch};
Dsave(1) = {Delta};
Rsave(1) = {Report};

% initial data - inner layer
adof = find(U);
nvox = numel(adof);
usum = [round(VolND*full(Notch(adof,1))'); ...
        round(VolND*full(Delta(adof,1))'); ...
        round(VolR*full(Report(adof,1))')];
umod.u0 = zeros(7,Nvoxels_,nvox);
for i = 1:nvox
  ind = floor(Nvoxels_*rand(1,usum(1,i)))+1;
  umod.u0(1,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
  ind = floor(Nvoxels_*rand(1,usum(2,i)))+1;
  umod.u0(2,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
  ind = floor(Nvoxels_*rand(1,usum(3,i)))+1;
  umod.u0(3,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
end

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

    % determine a suitable micro time-step
    if k > numel(Etime) || Etime(k) > tspan(j)
      ttend = tspan(j);
    else
      ttend = Etime(k);
    end

    % determine a suitable micro time-step
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

    % take K micro-steps
    for kk = 1:K
      % update junctional and protrusional singalling contacts
      umod.u0(4,:,:) = reshape(repmat(usum(1,:)*Nja'-usum(1,:),Nvoxels_,1),1,Nvoxels_,[]);
      umod.u0(5,:,:) = reshape(repmat(usum(1,:)*Npa'-usum(1,:),Nvoxels_,1),1,Nvoxels_,[]);
      umod.u0(6,:,:) = reshape(repmat(usum(2,:)*Nja'-usum(2,:),Nvoxels_,1),1,Nvoxels_,[]);
      umod.u0(7,:,:) = reshape(repmat(usum(2,:)*Npa'-usum(2,:),Nvoxels_,1),1,Nvoxels_,[]);
      % (in reality this belongs to the DLCM-layer: somewhat an abuse of the
      % RDME-layer)

      % parse/compile/run NSM-solver:
      umod.seed = umod.seed+1;
      umod = urdme(umod);
      % will not parse/compile again:
      umod.parse = 0;
      umod.compile = 0;
      
      % prepare for next step
      u = reshape(umod.U(:,end,:),[],Nvoxels_,ncell);
      usum = permute(sum(u,2),[1 3 2]);
      umod.u0 = u;
      % (note: umod.u0(4:7,:) is next overwritten)
    end

    % update state
    Notch(adof1,1) = sparse(usum(1,1:nvox)/VolND)';
    Delta(adof1,1) = sparse(usum(2,1:nvox)/VolND)';
    Report(adof1,1) = sparse(usum(3,1:nvox)/VolR)';
    Notch(adof2,2) = sparse(usum(1,nvox+1:end)/VolND)';
    Delta(adof2,2) = sparse(usum(2,nvox+1:end)/VolND)';
    Report(adof2,2) = sparse(usum(3,nvox+1:end)/VolR)';
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
      ID(n,1) = ID(ix,1); ID(ix,1) = 0;
      % change inner state in the same way
      adof_ = nonzeros(ID);
      u = u(:,:,adof_);
      usum = permute(sum(u,2),[1 3 2]);
      umod.u0 = u;
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
        ID(n,1) = ID(ix,m);
      else
        % one cell moved into an already occupied voxel
        Delta(n,2) = Delta(ix,m);
        Notch(n,2) = Notch(ix,m);
        Report(n,2) = Report(ix,m);
        ID(n,2) = ID(ix,m);
      end
      if m == 1
        Delta(ix,1) = Delta(ix,2);
        Notch(ix,1) = Notch(ix,2);
        Report(ix,1) = Report(ix,2);
        ID(ix,1) = ID(ix,2);
      end
      Delta(ix,2) = 0;
      Notch(ix,2) = 0;
      Report(ix,2) = 0;
      ID(ix,2) = 0;
      % change inner state in the same way
      adof_ = nonzeros(ID);
      u = u(:,:,adof_);
      usum = permute(sum(u,2),[1 3 2]);
      umod.u0 = u;
    elseif Etype(k) == 3
      % proliferation event
      ix = find(Event(:,k) == 1);
      % share Delta-Notch-Reporter between the two cells
      r1 = rand;
      Delta(ix,2) = r1*Delta(ix,1); Delta(ix,1) = (1-r1)*Delta(ix,1);
      r2 = rand;
      Notch(ix,2) = r2*Notch(ix,1); Notch(ix,1) = (1-r2)*Notch(ix,1);
      r3 = rand;
      Report(ix,2) = r3*Report(ix,1); Report(ix,1) = (1-r3)*Report(ix,1);
      ID(ix,2) = ID(ix,1);
      % change inner state in the same way
      adof_ = nonzeros(ID);
      u = u(:,:,adof_); % increases the cell population
      ncell = ncell+1;
      ix_ = find(adof_ == ID(ix,1));
      % share inner state analogously
      u(1:3,:,ix_(1)) = round(tprod(u(1:3,:,ix_(1)), ...
                                    [r1 r2 r3],[1 2],[3 1]));
      u(1:3,:,ix_(2)) = u(1:3,:,ix_(2))-u(1:3,:,ix_(1));
      usum = permute(sum(u,2),[1 3 2]);
      umod.u0 = u;
    end
    % finally, here is the actual event:
    U = U+Event(:,k);
    ID = spreplace(ID,1:ncell);
    k = k+1;
  end
  
  % record values
  Usave(j) = {U};
  Nsave(j) = {Notch};
  Dsave(j) = {Delta};
  Rsave(j) = {Report};
end

% display final frame
figure(1), clf,
i = numel(Usave);
NDRplot(Dsave{i},Usave{i},V,R);

return;

% movie visualization
mx = full(max(Dsave{end}(:)));
figure(1),
M = struct('cdata',{},'colormap',{});
for i = 1:numel(Usave)
  clf,
  NDRplot(Dsave{i},Usave{i},V,R,mx);
  M(i) = getframe(gcf);
end

movie2gif(M,{M([1:2 end]).cdata}, ...
          'visual/NDR_RDME_dynamic.gif', ...
          'delaytime',0.1,'loopcount',0);

return;
