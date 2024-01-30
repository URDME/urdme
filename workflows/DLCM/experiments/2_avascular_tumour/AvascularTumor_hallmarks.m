% Simulation of an avascular tumour model.
%
%   Avascular tumour growth: An initial circular population cells (one
%   per voxel) lie in a domain rich in oxygen. Cells consume oxygen at
%   a constant rate, lambda. Cells occupying a voxel with oxygen above
%   cutoff_prol can proliferate at a rate r_prol. Cells occupying
%   voxels with an oxygen concentration below cutoff_die can die at a
%   rate r_die.  Dead cells are represented with a voxel with value
%   -1, these dead cells can degrade and stop occupying space at a
%   rate r_degrade. Default setting disregards effects from oxygen, but
%   nonzero kappa_prol and kappa_death sets oxygen thresholds for
%   cell proliferation and death.
%
%   No oxygen consumption by dying cells. Dying cells constitute
%   pressure sinks due to mass loss. There is surface tension on 
%   every contour that is large enough.
%
%   Includes a minimal interpretation of the first four hallmarks
%   of cancer (and their absence) with local mutations at the cellular 
%   scale. Also includes a behaviour 'sluice'for each hallmark mutation
%   to activate first after a time hmk_dt*hmk#.
%
%   Permeability: Drate1 describes the rate diffusion rate of tumour
%   cells invading previously unvisited voxels. Drate2 is the rate
%   cells move into previously occupied but currently empty
%   voxels. Drate3 is the rate cells move into voxels that are already
%   occupied.

% E. Blom 2023-10-26 (hallmarks)
% E. Blom 2023-03-20 (revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

%% Setup & Simulation

% warnings
% noisy contour modified by Matlab at runtime, outcome not notably altered
warning('off', 'MATLAB:polyshape:repairedBySimplify')

% simulation interval
Tend = 300;
tspan_len = 91;
tspan = linspace(0,Tend,tspan_len);
report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

% The user specified cutoff and rate parameters for the proliferation,  
% death, degradation and consumption rules.
cons = 1;        % consumption of oxygen by cells
kappa_prol = 0;  % the minimum amount of oxygen for proliferation
mu_prol = 1;     % rate of proliferation of singly occupied voxels
kappa_death = 0; % the maximum amount of oxygen where cells can die
mu_death = 1;    % rate of death
mu_degrade = 0.1*mu_death; % rate of degradation for already dead cells
sigma = 1e-4;    % surface tension strength

% Savitzky-Golay filter parameters for contour-curvature evaluation
sgf_wind = 21;   % frame length
sgf_deg = 7;     % polynomial degree

% curvature hyperparameters
curv_cutoff = 0.1; % curvature larger than 1/curv_cutoff is clipped.
curv_thresh = 0.0; % minimum relative size of tumor clusters/holes

% Permeability parameters.
Drate1 = 1;     % into free matrix
Drate2 = 1;     % into already visited matrix
Drate3 = 25;    % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];

% hallmark parameters
min_gs_prol = 0.0001; % minimum growth signalling for cell prolif. (H1)
gs_decay = 500;       % growth signal decay (H1)
max_e_prol = 1;       % maximum occupied edge ratio for cell prolif. (H2)
mu_eps = 0.005;       % constant apoptosis rate factor for all cells (H3a)
mu_apoptosis = 0.25;  % conditioned apoptosis rate for all cells (H3b)
% baseline parameters:
hmk_baseline = [min_gs_prol, max_e_prol, mu_apoptosis]; % (H3c)
max_nprol = 20;       % initial number of proliferations per cell line (H4)

% hallmark mutation parameters
shift_global_hmk = true;   % gradually activate hmks 1-3 one-by-one
hmk_dt = 100;              % time period for each hmk activation
burnin = 0;                % burn-in time before activating hmks

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 51; % odd so the BC for oxygen can by centered

% all parameters
allpar = [cons kappa_prol mu_prol kappa_death mu_death mu_degrade ...
    sigma Drate_(:)'];

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');
[Ne, ~] = dt_neighe(V,R); % edge length operator
% Getting h for square discr. only
hmax = 2/(Nvoxels-1);

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[~,dM,N] = dt_operators(P,T);
[L,M] = assema(P,T,1,1,0); % get unscaled L
neigh = full(sum(N,2));

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));
[extdof_neigh, ~] = find(N(:, extdof)); % dofs neighbouring extdofs

% Initial population
IC = 1;
R0 = 1;   % initial radius for IC 1
R1 = 0;
R2 = 0;
U = set_init(IC,R1,R2,P,Nvoxels,R0);
U(extdof) = 0;  % avoid boundary issues
% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% assign unique cell indices
Uidx = fsparse( zeros(numel(U), 2)); 
u1 = find(U);                   % voxels with one cell
u2 = find(U==2);                % voxels with 2 cells
Uidx(u1, 1) = 1:numel(u1);      % assign indices for first cells in voxels
Uidx(u2, 2) = numel(u1) + (1:numel(u2));  % ... and second one, if exists
next_Uidx = numel(u1) + numel(u2) + 1;    % ensure next cell has unique idx
unused_Uidx = [];                         % reuse indices from dead voxels

% assign cell phenotypes
Utype = fsparse(zeros(4, numel(U)*2));
Utype(1, 1:next_Uidx-1) = max_nprol;        % nr of cell divisions left
Utype(2, 1:next_Uidx-1)  = min_gs_prol;     % growth signal req. for prol.
Utype(3, 1:next_Uidx-1)  = max_e_prol;      % occupied neighbour inhibition
Utype(4, 1:next_Uidx-1)  = mu_apoptosis;    % general apoptosis rate

% local hallmark mutation parameters
% row corresponds to Utype row phenotype
mut_rate = [0; [0.0; 0.0; 0.0]]; % mutation rate used from start
alpha_mut = [0.2 0.2 0.2];       % hallmark 1-3 mutation rate value
tel_mut_rate = 1e-3;  % mutation rate for large increase in n_prol
eta_tel = 1000;       % proliferation number increase by telomere mutation

ncells = sum(U(U>=0));  % debugging counter

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;
Prsave = cell(1,numel(tspan));   % save pressure distribution
Prsave{1} = 0.*U;                % save inital pressure state as zero
Oxysave = cell(1,numel(tspan));  % save pressure distribution
Oxysave{1} = 0.*U + 1;           % save inital ozygen state as one
Adofsave = cell(1,numel(tspan)); % save Adofs for Prsave plotting
Adofsave{1} = 1:length(U);
BFTsave = cell(1,numel(tspan));  % Fourier Transform of boundary contour
hmksave = cell(1,numel(tspan));  % save hallmark mean state per voxel

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

% growth signal operator
M = M - OLai*M + OLai; % Dirichlet at extdofs
GLa.X = (gs_decay*M + OLa.X);
[GLa.L,GLa.U,GLa.p,GLa.q,GLa.R] = lu(GLa.X,'vector');

% initialise reproduction number arrays
R_mean = [];
R_std = [];
R_dt = [];
 
ndof = int16.empty;
ndof_ = int16.empty; 
ddof = (sqrt(P(1,:).^2 + P(2,:).^2) < 0.2)';    % death dofs
Rdof = (sqrt(P(1,:).^2 + P(2,:).^2) < 0.6)';    % dofs used to measure R
% event counter
e_count = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);
Nesave = cell(1,numel(tspan));  % nr of events per type (birth, death, ...)
Nesave(1) = {e_count};
while tt <= tspan(end)
  % classify the DOFs
  adof = find(U); % all filled voxels
  % singularly occupied voxels on the boundary:
  bdof_m_both = find(N*(U ~= 0) < neigh & abs(U) == 1);
  bdof_m_live = find(N*(U ~= 0) < neigh & U == 1);
  % singly occ. on boundary to extdof are immovable - inf. tissue approx.
  bdof_m_both = setdiff(bdof_m_both, extdof_neigh); 
  bdof_m_live = setdiff(bdof_m_live, extdof_neigh);
  sdof = find(U > 1); % voxels with 2 cells
  ldof = find(U > 0); % voxels with living cells
  
  % voxels with 2 cells in them _which may move_, with a voxel
  % containing less number of cells next to it (actually 1 or 0):
  sdof_m = find(N*(U > 1) < neigh & U > 1);
  Idof = (N*(U ~= 0) > 0 & U == 0); % empty voxels touching occupied ones
  idof1 = find(Idof & ~VU); % "external" OBC1
  idof2 = find(Idof & VU);  % "internal" OBC2
  idof = find(Idof);

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];

  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';  
  [bdof_m_both_,bdof_m_live_,sdof_,sdof_m_,ldof_,idof1_,idof2_,idof_,... 
      adof_] = map(Adof_,Adof,bdof_m_both,bdof_m_live,sdof,sdof_m,...
      ldof,idof1,idof2,idof,adof);
    
  % update Laplacian if tumor domain changed
  if updLU
    % get curvature for large enough contours
    [X,Y,Z] = meshdata(P(1,:),P(2,:),abs(U)'); % U < 0 thus included 
    [cpos, ~] = contour(X, Y, Z, [0.5, 0.5]);    

    % get all contour segments
    boundaries = cell(1,0);
    cn = [];
    idx = 1;   % starting at index 1
    while idx < length(cpos)
      cn = [cn cpos(2,idx)];
      % save contour points
      boundaries{end+1} = [cpos(1,idx+1:idx+cn(end)); ...
                          cpos(2,idx+1:idx+cn(end))];
      idx = idx + cn(end) + 1;    % step to next contour line
    end

    % get curvature at idofs that are not alone
    desired_coords = [P(1,idof); P(2,idof)]; 
    curvature = curvature2D(boundaries, desired_coords, curv_cutoff, ...
        curv_thresh, sgf_wind, sgf_deg);

    % get main boundary for evaluating small perturbations
    cidx = find(cn == max(cn));
    boundary_main = boundaries{cidx};  % main contour line

    % set Dirichlet only at idofs that are from large enough contours
    idof_keep = find(~isnan(curvature));
    curvature = curvature(idof_keep);

    % pressure Laplacian
    La.X = L(Adof,Adof);
    % rescale 'true' idofs for Dirichlet BC
    Lai = fsparse(idof_(idof_keep),idof_(idof_keep),1,size(La.X));

    % however, rescale all idof hats for Neumann condition there
    Lai2 = fsparse(idof_,idof_,1,size(La.X));  
    % Scale Laplacian on boundary
    La.X = La.X - diag(sum(Lai2*La.X,2)); % replace with scaled hats (all)
    La.X = La.X-Lai*La.X+Lai; % remove fully supported hat functions

    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

    updLU = false; % assume we can reuse
  end
  
  % Oxygen on global domain
  Oxy = full(fsparse([extdof; adof],1, ...
                     [ones(size(extdof)); ... 
                      -cons*full(max(U(adof),0).*dM(adof))], ...
                     [size(OLa.X,1) 1]));
  Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));
  
  % find proliferating cells & dying cells (for pressure distr.)
  pdof = find(U > 0 & (Oxy > kappa_prol));     % prolif. dofs
  % do not use pdof for oxygen consumption!
  ndof = find(U==-1); 
  
  % extdofs not counted as idofs
  notidof = intersect(extdof, idof(idof_keep));
  
  % Do mapping again
  [ndof_,pdof_,notidof_] = ...          
  map(Adof_,Adof,ndof,pdof,notidof);

  % Pressure on local domain
  Pr = full(fsparse(sdof_,1, ...
                    [mu_prol.*dM(sdof)], ...
                    [size(La.X,1) 1]));     % RHS first...
  Pr(ndof_) = -mu_death.*dM(ndof);          % sinks
  Pr(idof_(idof_keep)) = +sigma.*curvature; % set BC
  Pr(notidof_) = 0;   %no surface tension on extdof - inf. tissue approx.
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution
  
  % Proliferation signal on global domain
  idof3 = find(VU & U == 0);    % all visited holes are sources
  Growth_signal = full(fsparse(idof3,1, ...
                    [dM(idof3)], ...
                    [size(GLa.X,1) 1]));     % sources are visited holes
  Growth_signal(extdof) = 0;                 % zero 'far-away'
  Growth_signal(GLa.q) = GLa.U\(GLa.L\(GLa.R(:,GLa.p)\Growth_signal));
  
  % intensities of possible events
  
  % (1a) moving boundary DOFs, [1, -1]->0
  [ii,jj_] = find(N(bdof_m_both,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) == 0);       % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(bdof_m_both_(ii))-Pr(jj_),0).* ...
                 Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                 numel(bdof_m_both));
  % rates scaled like FV method for cartesian grid!
  moveb0 = full(gradquotient*grad./hmax^2);
  % can generalize to scale by voxel volume here?
  
  % (1b) moving boundary DOFs, 1->1 (no necrotic onto single cell)
  
  % first find live bdofs that are only on tumor boundary
  [bdof_ii, ~] = find(N(bdof_m_live,idof(idof_keep)));  
  bdof_m_live = bdof_m_live(unique(bdof_ii));
  
  % Do mapping again
  moveb1 = double.empty([0,1]); % set to zero if no such dofs exist
  if ~isempty(bdof_m_live)
      [bdof_m_live_] = ...          
      map(Adof_,Adof,bdof_m_live);

      [ii,jj_] = find(N(bdof_m_live,Adof)); % neighbours...
      keep = find(U(Adof(jj_)) == 1);       % ...to move to
      ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);

      % remove any possibly remaining negative rates
      grad = fsparse(ii,1,max(Pr(bdof_m_live_(ii))-Pr(jj_),0).* ...
                     Drate_(1), ... % ensure consistent tumor front speed
                     numel(bdof_m_live));
                 % 1->1 for bdofs move at same Drate as 1->0 (non-visited)
      moveb1 = full(gradquotient*grad./hmax^2);
  end
  
  % (2) also certain sources may move by the same physics
  [ii,jj_] = find(N(sdof_m,Adof));                     % neighbours...
  keep = find(-1 < U(Adof(jj_)) & U(Adof(jj_)) < 2);   % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  jj0_ = VU(Adof(jj_)); % jj_ visited

  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0).* ...
               Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1), ...
                 numel(sdof_m)); % (abs as U could be -1)
  moves = full(gradquotient*grad./hmax^2); 
  
  % (3) proliferation/death/degradation rates
  H1 = fsparse([]);
  H2 = fsparse([]);
  H3 = fsparse([]);
  H4 = fsparse([]);

  H1(Adof_,1) = 0;                     % hallmark conditions that...
  H1(adof_,1) = Utype(2,Uidx(adof,1)); % properly map cell type to Adof
  
  H2(Adof_,1) = 0;
  H2(adof_,1) = Utype(3,Uidx(adof,1));
  
  H3(Adof_,1) = 0; 
  H3(adof_,1) = Utype(4,Uidx(adof,1));
  
  % H1, H4 here only care about singly occ.
  H4(Adof_,1) = 0;
  H4(adof_,1) = Utype(1,Uidx(adof,1));
  
  % extdofs are considered occupied for edge evaluation
  U_tmp = U; U_tmp(extdof) = 1;
  e_cell = Ne(Adof,Adof)*(U_tmp(Adof).*(U_tmp(Adof)>=0));
     
  birth = full(mu_prol*(U(Adof) == 1).*(Oxy(Adof) > kappa_prol)...
        .*(H4>0) ...                                            % H4
        .*(e_cell < H2(Adof_,1))...                             % H2
        .*(Growth_signal(Adof)>=H1(Adof_,1)));                  % H1

  total_birth = sum(birth); 
  birth(isnan(birth)) = 0;
  % (as we get some 0/0 terms if total_birth == 0);
 
  % hallmark conditions for 2nd cell in each voxel
  H1(Adof_,2) = 0;
  H1(sdof_,2) = Utype(2,Uidx(sdof,2));
  
  H2(Adof_,2) = 0; 
  H2(sdof_,2) = Utype(3,Uidx(sdof,2));
  
  H3(Adof_,2) = 0;
  H3(sdof_,2) = Utype(4,Uidx(sdof,2));  
  
  H4(Adof_,2) = 0;
  H4(sdof_,2) = Utype(1,Uidx(sdof,2)); 
         
  death2 = full((U(Adof)>0)... % hypoxia & imposed death at ddof
            .*((mu_death.*(Oxy(Adof) < kappa_death | ddof(Adof))) ... 
            + H3(Adof_,:).*( mu_eps ...        % H3a (uniform apoptosis)
            + (U(Adof) > 1).*((Growth_signal(Adof) < H1(Adof_,:)) ... % H3b
            + (e_cell >= H2(Adof_,:)) )... 
            ... % H3c, DNA damage per cell
            + (sqrt((1 - H1(Adof_,:)./hmk_baseline(1)).^2 ...
            + (1 - H2(Adof_,:)./hmk_baseline(2)).^2)) ))...
            ... % disregard non-doubly occ. in 2nd column
            .*[U(Adof)>0 U(Adof)>1] ); 
  death = sum(death2,2); % sum over doubly occupied voxels
  degrade = full(mu_degrade*(U(Adof) == -1));

  intens = [moveb0; moveb1; moves; birth; death; degrade];
  lambda = sum(intens);
  dt = -reallog(rand)/lambda; 
  rnd = rand*lambda;
  cum = intens(1);
  ix_ = 1;
  while rnd > cum
    ix_ = ix_+1;
    cum = cum+intens(ix_);
  end
  % (now ix_ points to the intensity which fired first)

  % Metric from branching process
  ldof2 = intersect(ldof,find(Rdof));   % living dofs in metric subset
  sdof2 = intersect(sdof,find(Rdof));
  [ldof2_] = ...          
  map(Adof_,Adof,ldof2);
  Usub = U(Rdof);                       % cells in metric subset
  Usub_tot = sum(Usub(Usub>0)); 
  Q = birth(ldof2_)./(death(ldof2_)+birth(ldof2_)); % process 'success'
  R_i = (Q-Q.^(Utype(1,Uidx(ldof2,1))'+1))./(1-Q);  % cell reproduction nr.
  % sdofs have metric = 0 so divide mean and var by Usub_tot
  R_mean = [R_mean 1/Usub_tot*sum(R_i)];     
  R_std = [R_std sqrt((sum((R_i-R_mean(end)).^2) ...
      + 2*numel(sdof2)*R_mean(end)^2)/Usub_tot)];
  R_dt = [R_dt tt];
   
  % report back
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    Prsave(i+1:iend) = {Pr};
    Oxysave(i+1:iend) = {Oxy};
    Adofsave(i+1:iend) = {Adof};
    Nesave(i+1:iend) = {e_count};
    
    % hallmark mean at each voxel
    hmk_mean = (abs(1 - H1(Adof_,:)./hmk_baseline(1)) ...
    + abs(H2(Adof_,:)./hmk_baseline(2)-1) ...
    + abs(1 - H3(Adof_,:)./hmk_baseline(3)) ...
    + H4(Adof_, :)./(H4(Adof_,:) + max_nprol)).*[U(Adof)>0 U(Adof)>1];
    % take cell mean over each doubly occ. voxel
    hmk_mean(hmk_mean(:,2) > 0,1) = mean(hmk_mean(hmk_mean(:,2) > 0, :),2);
    hmk_mean = hmk_mean./4; % div. by nr of hallmarks
    hmksave(i+1:iend) = {hmk_mean(:,1)};
    
    % population mean of hallmark variables
    gs_prol(i+1:iend) = mean([Utype(2,Uidx(U>0,1)) Utype(2,Uidx(U>1,2))]);
    e_prol(i+1:iend) = mean([Utype(3,Uidx(U>0,1)) Utype(3,Uidx(U>1,2))]); 
    mu_apt_all(i+1:iend) = mean([Utype(4,Uidx(U>0,1)) ... 
        Utype(4,Uidx(U>1,2))]);
    nprol_avg(i+1:iend) = ...
    mean([Utype(1,Uidx(U>0,1)) Utype(1,Uidx(U>1,2))]);

    % monitor the maximum outlier cell:
    max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));

    % the number of occupied voxels
    num_vox = sum(VU);
    
    % 'radius' of total volume, pi*r^2=V
    volrad = sqrt(num_vox*hmax^2/pi);

    % the rates
    inspect_rates = [sum(moveb0) sum(moveb1) sum(moves) ...
                     sum(birth) sum(death) sum(degrade)];
                 
    cr = sqrt(boundary_main(1,:).^2 + boundary_main(2,:).^2);

    % estimate contour roundness from its corresponding polygon
    pgon = polyshape(boundary_main(1,:),  boundary_main(2,:));
    perim = perimeter(pgon);
    parea = area(pgon);
    roundness(i+1:iend) = perim^2/(4*pi*parea);  % = 1 for perfect circles

    % save the Fourier transform of the contour
    Yfft = fft(cr);
    Fs = numel(cr);
    n = Fs;
    P2 = abs(Yfft/n);
    n2 = floor(n/2); % avoid fractional indexing
    P1 = P2(:, 1:n2 + 1);
    P1(:, 2:end-1) = 2*P1(:, 2:end-1);
    BFTsave(i+1:iend) = {[0:1:(n2-1);P1(1:n2)]};
                 
    % P:Q:N ratios (areas reduced to radii)
    % tumor radius
    rp(i) = sqrt(numel(adof)*hmax*hmax/pi);
    % quiescent-prolif. interface radius
    rq(i) = sqrt((numel(adof) - numel(pdof))*hmax*hmax/pi);     % wrt pdof
    % necrotic-quiescent interface radius
    dying_dof = find(U & (Oxy < kappa_death));
    rn(i) = sqrt(numel(ndof)*hmax*hmax/pi);
    
    i = iend;
    
  end

  if ix_ <= numel(moveb0) % [1, -1] -> 0
    e_count.moveb = e_count.moveb+1;
    % movement of a boundary (singly occupied) voxel
    ix_ = bdof_m_both_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    % (will only move into an empty voxel:)
    jx_ = jx_(U(Adof(jx_)) == 0);
    rates = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    U(n) = U(ix);             % note: U(ix) can be both -1 and 1 here!
    U(ix) = 0;
    Uidx(n,1) = Uidx(ix, 1);  % move cell index
    Uidx(ix, 1) = 0;
    updLU = true;             % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1) % 1 -> 1
    e_count.moveb = e_count.moveb+1;
    % movement of a boundary (singly occupied) voxel
    ix_ = ix_-numel(moveb0);
    ix_ = bdof_m_live_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    % (will only move into singly occupied voxel:)
    jx_ = jx_(U(Adof(jx_)) == 1);
    rates = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    U(n) = U(n)+1;
    U(ix) = U(ix)-1;
    Uidx(n, 2) = Uidx(ix, 1);    % move index to occupied cell
    Uidx(ix, 1) = 0;
    if U(ix) < 0
        error('movement created necrotic cell!')
    end
    updLU = true; % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves) % 2 -> [0,1]
    e_count.moves = e_count.moves+1;
    % movement of a cell in a doubly occupied voxel
    ix_ = ix_-numel(moveb0)-numel(moveb1);
    ix_ = sdof_m_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    % (won't move into a voxel containing a dead -1 cell:)
    jx_ = jx_(-1 < U(Adof(jx_)) & U(Adof(jx_)) < 2);
    rates = Drate_(2*VU(Adof(jx_))+abs(U(Adof(jx_)))+1).* ...
            max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    if U(n) == 0, updLU = true; end % boundary has changed
    U(n) = U(n)+1;
    U(ix) = U(ix)-1;
    Uidx(n,U(n)) = Uidx(ix, 1);  % first in - last out
    Uidx(ix, 1) = Uidx(ix, 2);
    Uidx(ix, 2) = 0;
    if U(ix) < 0
        error('movement created necrotic cell!')
    end
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves)+numel(birth)
    e_count.birth = e_count.birth+1;
    % proliferation
    birth_count = birth_count+1;
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves);
    ix = Adof(ix_);
    U(ix) = U(ix)+1;
    % assign index to new cell
    if isempty(unused_Uidx)
        Uidx(ix, 2) = next_Uidx;
        next_Uidx = next_Uidx + 1;
    else
        Uidx(ix, 2) = unused_Uidx(end);     
        unused_Uidx = unused_Uidx(1:end-1);    % index is now used
    end
    % update (or mutate) current cell phenotype
    Utype(1, Uidx(ix,1)) = Utype(1, Uidx(ix,1)) - 1;   % telomere reduction
    dUtype = -ones(numel(mut_rate),2);
    while sum(dUtype(:) < 0) || sum(dUtype(3,:) > 2.01)% reject bad values
        dUtype(:, 1) = Utype(:, Uidx(ix,1)).*(1 ...
            + mut_rate.*randn(numel(mut_rate),1));
        dUtype(:, 2) = Utype(:, Uidx(ix,1)).*(1 ...
            + mut_rate.*randn(numel(mut_rate),1));
        dUtype(3, :) = Utype(3, Uidx(ix,1)).*(1 ...
            + mut_rate(3).*randn(1,2)...
            .*(2.01-Utype(3, Uidx(ix,1))));
    end
    Utype(:, Uidx(ix,1)) = dUtype(:, 1);
    Utype(:, Uidx(ix,2)) = dUtype(:, 2);
    Utype(1, Uidx(ix,:)) = Utype(1, Uidx(ix,:)) ...
        + eta_tel*(rand(1,2) < tel_mut_rate);
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves)+numel(birth)...
          +numel(death)
    e_count.death = e_count.death+1;
    % death
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves)-numel(birth);
    ix = Adof(ix_);
    if U(ix) == 2
      U(ix) = 1; % (removed directly)
      % find which cell died and update indices
      rates = mu_death.*(Oxy(ix) < kappa_death)... 
            + H3(ix_,:) .* (mu_eps ...  
            + ((Growth_signal(ix) < H1(ix_,:)) ...
            + (e_cell(ix_) >= H2(ix_,:))) ...
            + (sqrt((1 - H1(ix_,:)./hmk_baseline(1)).^2 ...
            + (1 - H2(ix_,:)./hmk_baseline(2)).^2)) );
      m = find(cumsum(rates) > rand*sum(rates),1,'first');
      unused_Uidx = [unused_Uidx Uidx(ix, m)];    % free index for reuse
      Uidx(ix, 1) = Uidx(ix, mod(m,2)+1);
      Uidx(ix, 2) = 0;
      e_count.degrade = e_count.degrade+1;
    else
      U(ix) = -1;
    end
  else
    e_count.degrade = e_count.degrade+1;
    % degradation
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves)-numel(birth)...
        -numel(death);
    ix = Adof(ix_);
    U(ix) = 0;
    unused_Uidx = [unused_Uidx Uidx(ix, 1)];    % free index for reuse
    Uidx(ix, 1) = 0;
    updLU = true; % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  end
  U(extdof) = 0; % avoid boundary issues for spatially 'infinite' tissue
  % free index for reuse
  unused_Uidx = [unused_Uidx, (Uidx(extdof,1).*(Uidx(extdof,1)>0))'];    
  unused_Uidx = unused_Uidx(unused_Uidx ~= 0);
  Uidx(extdof, 1) = 0;
  
  % update metric for cell number consistency
  ncells = sum(U(U>=0));
  
  if shift_global_hmk
  % brute-force initiate hallmarks
  if tt > burnin   % wait until after burn-in period
    mut_rate(2) = alpha_mut(1);     % allow HMK1 mutation
  end
  if tt > burnin + hmk_dt
      mut_rate(3) = alpha_mut(2);   % allow HMK2 mutation
  end
  if tt > burnin + 2*hmk_dt
      mut_rate(4) = alpha_mut(3);   % allow HMK3 mutation
  end
  end
  
  tt = tt+dt;
  report(tt,U,'');
    
  % update the visited sites
  VU = VU | U;
end
report(tt,U,'done');

%% Plots

% plot images at times tstamps
if numel(tspan) > 80
 tstamps = [4 37 46 64 70 79];  % default interest points
else
 tstamps = 2:round(numel(tspan)/6):numel(tspan); % evenly distributed
end

% get VU for each tstamp
VUsave = cell(1,numel(tspan));
VU_tmp = (Usave{1} ~= 0);
for i = 2:numel(tspan)
    VU_tmp = VU_tmp | Usave{i};
    VUsave{i} = VU_tmp;
end

M = struct('cdata',{},'colormap',{});
figure(3), clf,
t = tiledlayout(2,3, 'TileSpacing','none','Padding','none');
im = 1;

i = 1;
for iii = tstamps

  outdof = find(VUsave{iii} & Usave{iii} == 0);% previously visited voxels
  
  % finding Adofs from Usave
  adof = find(Usave{iii}); % all filled voxels
  % empty voxels touching occupied ones
  Idof = (N*(Usave{iii} ~= 0) > 0 & Usave{iii} == 0); 
  idof1 = find(Idof & ~VUsave{iii}); % "external" OBC1
  idof2 = find(Idof & VUsave{iii});  % "internal" OBC2
  idof = find(Idof);
 % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];
  
  figure(3)
  nexttile
  c = [0.6 0.6 0.6];   % darker grey than background
  axis(1.*[-1 1 -1 1]); axis square, axis off
  patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');
  patch('Faces',R(idof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');
  
  patch('Faces',R(Adof,:),'Vertices',V,'FaceVertexCData', ...
        full(hmksave{iii}),'FaceColor','flat', 'EdgeColor', 'None');
  patch('Faces',R(idof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');
  try
    clim([0,1]);
  catch
    caxis([0,1]);
  end
  drawnow;
  
  % uncomment to save gif (1/2)
  %M(im) = getframe(gcf);
  %im = im+1;
  
  i = i+1;
end

figure(3)
set(gcf,'Position',[100 100 340 240]);
cb = colorbar;
cb.Layout.Tile = 'east';

% uncomment to save gif  (2/2)
%movie2gif(M,{M([1:2 end]).cdata},'tumor_animation.gif', ...
%          'delaytime',0.1,'loopcount',0);

% R moving average + std
figure(6)
k = 20;
% get moving avg & std
R_mov = movmean(R_mean, k, 'SamplePoints', R_dt);
R_movstd = movstd(R_mean, k, 'SamplePoints', R_dt);
plot(R_dt(1:1000:end), R_mov(1:1000:end))
hold on
errorshade(R_dt(1:1000:end), R_mov(1:1000:end) - R_movstd(1:1000:end), ...
R_mov(1:1000:end) + R_movstd(1:1000:end))
ylabel('$R_\textrm{pop}$', 'Interpreter','latex')
axis([0, tspan(end), 0, 5])
hold on
% plot hallmarks
yyaxis right
plot(tspan(2:end),1-gs_prol(2:end)./gs_prol(2), 'Linewidth', 1)
plot(tspan(2:end),e_prol(2:end)./e_prol(2)-1.0, 'Linewidth', 1)
plot(tspan(2:end),1-mu_apt_all(2:end)./mu_apt_all(2), 'Linewidth', 1)
plot(tspan(2:end),(nprol_avg(2:end))./(nprol_avg(2:end)+nprol_avg(2)), ...
    'Linewidth', 1)
legend('Fitness', 'HMK1', 'HMK2', 'HMK3', 'HMK4')
grid on
xlabel('Time [$\mu_{\textrm{prol}}^{-1}$]', 'Interpreter','latex')
ylabel('', 'Interpreter','latex')
axis([0, tspan(end), -0.0, 1.2])
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 240]);
set(gca, 'fontname', 'Roman', 'FontSize', 9.0)
plot([tspan(tstamps); tspan(tstamps)],[0 5;0 5;0 5;0 5;0 5;0 5]', 'k--')
legend('Fitness', 'HMK1', 'HMK2', 'HMK3', 'HMK4')

return;
%% Save data & animate
 saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'Prsave', ...
     {Prsave}, 'Oxysave', {Oxysave}, 'Nesave', {Nesave}, ...
     'tspan', {tspan}, 'Adofsave', {Adofsave}, ...
     'R_dt', {R_dt}, 'R_mean', {R_mean}, 'R_std', {R_std}, ...
     'R', {R}, 'V', {V}, 'N', {N},  'hmksave', {hmksave}, ...
     'Utype', {Utype}, 'Uidx', {Uidx}, 'gs_prol', {gs_prol}, 'e_prol', ...
     {e_prol}, 'mu_apt_all', {mu_apt_all}, 'nprol_avg', {nprol_avg}, ...
     'Pr', {Pr}, 'Adof', {Adof},'adof', {adof}, 'adof_', {adof_}, ...
     'idof', {idof},'idof_', {idof_}, 'Nvoxels',{Nvoxels}, ...
     'P', {P}, 'E', {E}, 'T', {T}, 'roundness', {roundness}, 'radius', ...
     {R0}, 'rp', {rp}, 'rq', {rq}, 'rn', {rn}, 'hmax', {hmax}, ...
     'Tend', {Tend}, 'all_parameters', {allpar}, 'BFTsave', {BFTsave});
 filename = "AThallmarks_longrun";
 filename_saveData = filename + ".mat";
 save(filename_saveData,'-struct','saveData');
 
 return
 % animate tumor growth
 animate_growth(filename_saveData, "DLCM", true)
