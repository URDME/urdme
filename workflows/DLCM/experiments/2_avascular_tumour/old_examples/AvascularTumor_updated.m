% Simulation of an avascular tumour model.
%
%   Avascular tumour growth: An initial circular population cells (one
%   per voxel) lie in a domain rich in oxygen. Cells consume oxygen at
%   a constant rate, lambda. Cells occupying a voxel with oxygen above
%   cutoff_prol can proliferate at a rate r_prol. Cells occupying
%   voxels with an oxygen concentration below cutoff_die can die at a
%   rate r_die.  Dead cells are represented with a voxel with value
%   -1, these dead cells can degrade and stop occupying space at a
%   rate r_degrade.
%
%   No oxygen consumption by dying cells. Dying cells constitute
%   pressure sinks due to mass loss.
%   Surface tension on every contour that is 'large enough'
%
%   Permeability: Drate1 describes the rate diffusion rate of tumour
%   cells invading previously unvisited voxels. Drate2 is the rate
%   cells move into previously occupied but currently empty
%   voxels. Drate3 is the rate cells move into voxels that are already
%   occupied.

% E. Blom 2023-03-20 (revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

%%
% warnings
% noisy contour modified by Matlab at runtime, outcome not notably altered
warning('off', 'MATLAB:polyshape:repairedBySimplify')

% simulation interval
Tend = 20;
tspan_len = 101;
tspan = linspace(0,Tend,tspan_len);
report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

% The user specified cutoff and rate parameters for the proliferation,  
% death, degradation and consumption rules.
cons = 1;                  % consumption of oxygen by cells
kappa_prol = 0.94;         % minimum amount of oxygen for proliferation
mu_prol = 1;               % rate of proliferation of singly occupied voxels
kappa_death = 0.93;        % maximum amount of oxygen where cells can die
mu_death = 0.5;            % rate of death
mu_degrade = 0.1*mu_death; % rate of degradation for already dead cells
sigma = 0e-4;              % surface tension strength

% Savitzky-Golay Filter (sgf) parameters for curvature evaluation
sgf_wind = 21;      % frame length
sgf_deg = 7;        % polynomial degree

curv_cutoff = 0.1;  % curvature larger than 1/curv_cutoff is clipped.
curv_thresh = 0.95; % min. relative size of tumor clusters/holes
                    % w. surface tension

% Permeability parameters.
Drate1 = 1;         % into free matrix
Drate2 = 1;         % into already visited matrix 
Drate3 = 25;        % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 101; % odd so the BC for oxygen can by centered

% all parameters
allpar = [cons kappa_prol mu_prol kappa_death mu_death mu_degrade ...
    sigma Drate_(:)'];

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');
% Getting h for square discr. only
hmax = 2/(Nvoxels-1);

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[~,dM,N] = dt_operators(P,T);
[L,M] = assema(P,T,1,1,0);  % get unscaled L
neigh = full(sum(N,2));

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

% Initial population
IC = 1;
R0 = 0.1;   % initial radius for IC 1
R1 = 0;
R2 = 0;
U = set_init(IC,R1,R2,P,Nvoxels,R0);

ncells = sum(U(U>=0));  % debugging counter for mass imbalance checks

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

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

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);

% oxygen Laplacian - can construct this outside time loop!
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

ndof = int16.empty;
ndof_ = int16.empty;
% event counter
Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);
while tt <= tspan(end)
  % classify the DOFs
  adof = find(U); % all filled voxels
  % singularly occupied voxels on the boundary:
  bdof_m_both = find(N*(U ~= 0) < neigh & abs(U) == 1);
  bdof_m_live = find(N*(U ~= 0) < neigh & U == 1);
  sdof = find(U > 1); % voxels with 2 cells
  
  % voxels with 2 cells in them _which may move_, with a voxel
  % containing less number of cells next to it (actually 1 or 0):
  sdof_m = find(N*(U > 1) < neigh & U > 1);
  Idof = (N*(U ~= 0) > 0 & U == 0); % empty voxels touching occupied ones
  idof1 = find(Idof & ~VU);         % "external" OBC1
  idof2 = find(Idof & VU);          % "internal" OBC2
  idof = find(Idof);

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];

  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';  
  [bdof_m_both_,bdof_m_live_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_] ...
      = map(Adof_,Adof,bdof_m_both,bdof_m_live,sdof,sdof_m,idof1, ...
            idof2,idof,adof);
    
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
        idx = idx + cn(end) + 1;       % step to next contour line
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
    La.X = La.X-Lai*La.X+Lai;

    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

    updLU = false; % assume we can reuse
  end
  
  % RHS source term proportional to the over-occupancy and BCs
  Oxy = full(fsparse([extdof; adof],1, ...
                     [ones(size(extdof)); ... 
                      -cons*full(max(U(adof),0).*dM(adof))], ...
                     [size(OLa.X,1) 1]));
  Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));
  
  % find proliferating cells & dying cells (for pressure distr.)
  pdof = find(U > 0 & (Oxy > kappa_prol));     % prolif. dofs
  % do not use pdof for oxygen consumption!
  ndof = find(U==-1); 
  
  % Do mapping again
  [ndof_,pdof_] = ...          
  map(Adof_,Adof,ndof,pdof);

  % RHS source term proportional to the over-occupancy and BCs
  Pr = full(fsparse(sdof_,1, ...
                    [mu_prol.*dM(sdof)], ...
                    [size(La.X,1) 1]));     % RHS first...
  Pr(ndof_) = -mu_death.*dM(ndof);          % (sinks,
  Pr(idof_(idof_keep)) = +sigma.*curvature; % set BC)
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution
  
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
  birth = full(mu_prol*(U(Adof) == 1).*(Oxy(Adof) > kappa_prol));
  total_birth = sum(birth);
  birth(isnan(birth)) = 0;
  % (as we get some 0/0 terms if total_birth == 0);
 
  death = full(mu_death*(U(Adof) > 0).*(Oxy(Adof) < kappa_death));
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

  % report back
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    Prsave(i+1:iend) = {Pr};
    Oxysave(i+1:iend) = {Oxy};
    Adofsave(i+1:iend) = {Adof};

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
    roundness(i+1:iend) = perim^2/(4*pi*parea);  % = 1: perfect circle

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
    rq(i) = sqrt((numel(adof) - numel(pdof))*hmax*hmax/pi);
    % necrotic-quiescent interface radius
    dying_dof = find(U & (Oxy < kappa_death));
    rn(i) = sqrt(numel(ndof)*hmax*hmax/pi);

    i = iend;
  end

  if ix_ <= numel(moveb0) % [1, -1] -> 0
    Ne.moveb = Ne.moveb+1;
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
    U(n) = U(ix);   % note: U(ix) can be both -1 and 1 here!
    U(ix) = 0;
    updLU = true;   % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1) % 1 -> 1
    Ne.moveb = Ne.moveb+1;
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
    updLU = true; % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves)
    Ne.moves = Ne.moves+1;
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
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves)+numel(birth)
    Ne.birth = Ne.birth+1;
    % proliferation
    birth_count = birth_count+1;
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves);
    ix = Adof(ix_);
    U(ix) = U(ix)+1;
  elseif ix_ <= numel(moveb0)+numel(moveb1)+numel(moves)+numel(birth)+numel(death)
    Ne.death = Ne.death+1;
    % death
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves)-numel(birth);
    ix = Adof(ix_);
    if U(ix) == 2
      U(ix) = 1; % (removed directly)
      Ne.degrade = Ne.degrade+1;
    else
      U(ix) = -1;
    end
  else
    Ne.degrade = Ne.degrade+1;
    % degradation
    ix_ = ix_-numel(moveb0)-numel(moveb1)-numel(moves)-numel(birth)-numel(death);
    ix = Adof(ix_);
    U(ix) = 0;
    updLU = true; % boundary has changed
    if ncells ~= sum(U(U>=0))
        error('mass imbalance')
    end
  end
  ncells = sum(U(U>=0));
  tt = tt+dt;
  report(tt,U,'');
    
  % update the visited sites
  VU = VU | U;
end
report(tt,U,'done');

% when generating synthetic data, do not visualise
if exist('gendata_par', 'var')
    return;
end

%% Plots

% get VU for each tstamp if nonexistent
if  ~exist('VUsave', 'var')
    VUsave = cell(1,numel(tspan));
    VU_tmp = (Usave{1} ~= 0);
    for i = 2:numel(tspan)
        VU_tmp = VU_tmp | Usave{i};
        VUsave{i} = VU_tmp;
    end
end

if ~exist('cons', 'var')
    % load parameters
    cons = all_parameters(1);
    kappa_prol = all_parameters(2);
    mu_prol = all_parameters(3);
    kappa_death = all_parameters(4);
    mu_death = all_parameters(5);
    mu_degrade = all_parameters(6);
    sigma = all_parameters(7);
end

if ~exist('tstamps', 'var')
    % set default times for figures
    tstamps = [2 round(0.8*numel(tspan)) numel(tspan)];
end

i = 1; % figure index
for iii = tstamps
    fig=figure(i);
    set(gcf,'Position',[100 100 340/3 240/2]);
    clf,
  
    % prolif. region
    pdof = find(Usave{iii} > 0 & (Oxysave{iii} > kappa_prol)); 
    qdof = find(Usave{iii} > 0 & (Oxysave{iii} <= kappa_prol) ...
                           & (Oxysave{iii} > kappa_death)); % quiesc. region
    ndof = find(Usave{iii} == - 1);                         % necrotic cells
    % living, but dying
    ddof = find(Usave{iii} > 0 & (Oxysave{iii} <= kappa_death)); 
  
    outdof = find(VUsave{iii} & Usave{iii} == 0); % visited regions
  
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis(0.75.*[-1 1 -1 1]); axis square, axis off
  
    c = [0.6 0.6 0.6];   % darker grey than background
    patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
        'FaceColor','flat', 'EdgeColor','none');
 
    single = patch('Faces',R(qdof,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'), 'EdgeColor','none');
  
    double = patch('Faces',R(pdof,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'), 'EdgeColor','none');
  
    dead = patch('Faces',R(ndof,:),'Vertices',V, ...
        'FaceColor',[0 0 0], 'EdgeColor','none');
    
    dying = patch('Faces',R(ddof,:),'Vertices',V, ...
        'FaceColor',[0.25 0.25 0.25], 'EdgeColor','none');  
    
    % mark double
    ii = find(Usave{iii} == 2);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0], 'FaceAlpha', 0.45, 'EdgeColor','none');
    
    % mark domain origo by crosshairs
    plot([0 0], [-1 1],'w')
    plot([-1 1], [0 0] ,'w')
    
    drawnow;
  
    i = i+1;
end

% plot volumetric growth
PDE = false;    % plotting DLCM
RegionalCharacteristicsVisualiser

return;
%% Save data & animate
saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'Prsave', ...
                   {Prsave}, 'Oxysave', {Oxysave}, 'tspan', {tspan}, ...
                   'R', {R}, 'V', {V}, 'N', {N}, 'Pr', {Pr}, 'Adof', ...
                   {Adof},'adof', {adof}, 'adof_', {adof_}, 'idof', ...
                   {idof},'idof_', {idof_}, 'Nvoxels',{Nvoxels}, ...
                   'P', {P}, 'roundness', {roundness}, 'radius', ...
                   {R0}, 'rp', {rp}, 'rq', {rq}, 'rn', {rn}, 'hmax', ...
                   {hmax}, 'Tend', {Tend}, 'all_parameters', {allpar}, ...
                   'BFTsave', {BFTsave}, 'tstamps', {tstamps});
filename = "DLCMtumor_test";
filename_saveData = filename + ".mat";
save(filename_saveData,'-struct','saveData');
 
return
% animate tumor growth
animate_growth(filename_saveData, "DLCM", false)
