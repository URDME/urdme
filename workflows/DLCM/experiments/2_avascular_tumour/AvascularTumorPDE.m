% Simulation of a PDE avascular tumour model.
%
%   Avascular tumour growth: An initial density of cells lie in a domain
%   of oxygen. Cells consume oxygen at a constant rate, cons. 
%   Cells in a volume with oxygen above cutoff_prol proliferate at a rate
%   r_prol. Cells occupying volumes with an oxygen concentration below 
%   cutoff_die die at a rate r_die (and consume no oxygen). 
%   Mass balance together with 
%   cell flow velocity gives pressure sources and sinks corresponding to
%   growth and death respectively. 'Young-Laplace' surface tension
%   acts on on tumor boundary proportional to a constant sigma.
%
%   Solves oxygen and pressure equations by Finite Element Method.
%   Solves advection by Finite Volume method (cf. AvascularTumor.m).

% E. Blom 2023-01-01 (major: FV and curvature implementation)
% E. Blom 2022-08-01

% Dependencies
%  - RegionalDynamics2D
%  - solve_oxygen
%  - connect_boundary
%  - curvature2D
%  - dt_operators
%  - mesh2dual
%  - getmidpointcircle
  
%% Initial experiment setup

clear; clf;
 
% warnings
% noisy contour modified by Matlab at runtime, outcome not notably altered
warning('off', 'MATLAB:polyshape:repairedBySimplify')

seed = 314;            % for the gaussian noise
rng(seed)

% the user specified cutoff and rate parameters for the proliferation,
% death, degradation and consumption rules.
lambda = 1.15;         % consumption of oxygen by cells
kappa_prol = 0.94;     % the minimum amount of oxygen for proliferation
mu_prol = 1;           % rate of proliferation of singly occupied voxels
kappa_death = 0.93;    % the maximum amount of oxygen where cells can die
mu_death = 1.35;       % rate of death
sigma = 0e-4;          % surface tension coefficient
cutoff_bdof = 0.9;     % boundary density threshold
noise = 0.025;         % boundary flux noise factor
c_out = 1;             % oxygen source far away

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 101; % odd so the BC for oxygen can by centered

Tend = 10;                  % final time step

% Initial condition parameters
start_value = 1; % cell concentrations in the initial blob
radius = 0.1;
eps = 0.00;
k = 0;

% Curvature parameters
curv_thresh = 0.95;          % contour size cutoff
curv_cutoff = 2/(Nvoxels-1); % curvature > 1/curv_cutoff is clipped.

% Savitzky-Golay Filter (sgf) parameters for curvature evaluation
sgf_wind = 21;      % frame length
sgf_deg = 7;        % polynomial degree

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels); 
hmax = 2/(Nvoxels-1);
Nnodes = size(P,2);
[V,R] = mesh2dual(P,E,T, 'voronoi');

D = 1; % D_rate, the rate with which cells move in the domain. 

% simulation interval
n_measurements = 10*Tend+1;      % ~5-10 measurements per time unit is OK
tspan = linspace(0,Tend,n_measurements);
report(tspan,'timeleft','init'); %(this estimator gets seriously confused!)
dt = Tend/n_measurements;        % temporary dt

% store all parameters for when saving experiment results
allpar = [lambda, kappa_prol, mu_prol, kappa_death, mu_death, ...
          cutoff_bdof, sigma];

% initial population: circular blob of living cells, perturbed by mode k
r = sqrt( (P(1,:)).^2+(P(2,:)).^2);
theta = atan2(P(2,:),P(1,:));
ii = find(r < radius + eps.*cos(k*theta)); % radius of the initial blob

iinan = find(isnan(theta)); % find if any thetas are NaN...
ii = union(ii,iinan);       % ... assume they should be in IC also

%### comparison with radially symmetric analytical solution
% find stationary volume within given time frame
init_volume = radius^2*pi;
long_time = 100;
par = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
             'kappa_prol',kappa_prol,'kappa_death',kappa_death, ...
             'lambda',lambda,'c_out',c_out);
reg_char = RegionalDynamics2D(par, init_volume, linspace(0,long_time,1000)); 
Vn_max = max(reg_char.Vn);
Vq_max = max(reg_char.Vq);
Vp= reg_char.Vp;
Vp(end) = [];
Vp_max = max(Vp);
% evaluate stabilizing effect from sigma given expected stationary state
m = 1:20;  % check only up to mode 20
args = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
              'rp',sqrt(Vp_max/pi),'rq',sqrt(Vq_max/pi), ...
              'rn',sqrt(Vn_max/pi),'sigma',sigma,'modes',m) ;
sigma_req = get_growth_factors(args);

unstable_modes = sum(0 < sigma_req);   % sigma_req decreases with m
most_unstable = find(sigma_req == max(sigma_req));
disp('Note: with sigma = ' + string(sigma) + ...
     ', the tumor is unstable up to mode ' + string(unstable_modes) + ...
     ', and mode ' + string(most_unstable) + ' is most unstable');
%###

U = fsparse(ii(:),1,start_value, [Nnodes 1]); % intialize cell population

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[~,dM,~] = dt_operators(P,T);

[L,M] = assema(P,T,1,1,0); % assemble Laplacian and Mass matrix
N = (M ~= 0);              % Neighbouring matrix, with nonzero diagonal)

% find the external boundary dofs
% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

% connect boundary dofs
extdof = connect_boundary(extdof, M, N, P);

% visit marker matrix: 1 for voxels who have been occupied 
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;

% keeps track of the pressure for visuals
Prsave = cell(1,numel(tspan));   % save pressure distribution
Prsave{1} = 0.*U;                % save inital pressure state as zero
Adofsave = cell(1,numel(tspan)); % save Adofs for Prsave plotting
Adofsave{1} = 1:length(U);

% keeps track of oxygen and dofs for figures and understanding
Oxysave = cell(1,numel(tspan));

BFTsave = cell(1,numel(tspan)); % Fourier Transform of boundary contour

tt = tspan(1);
i = 1;

% LU-factorization structs
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
iOLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0); % for iterative oxygen solver

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   

% prepare oxygen 'pseudo-timestepping' machinery
M = M - OLai*M + OLai;
% rule of thumb: dtau < 1e-3
dtau = 2e-4;                  % time step for the local oxygen solver
iOLa.X = (M + dtau*OLa.X);    % inverse matrix for the local oxygen solver
[iOLa.L,iOLa.U,iOLa.p,iOLa.q,iOLa.R] = lu(iOLa.X,'vector');

[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

ndof = int16.empty; % no necrotic dofs in initial config. assumed!
pdof = find(U);     % all intial cells are proliferating, assumed
ndof_ = int16.empty;
pdof_ = (1:numel(pdof))'; 
% for solving pressure with 'far away' BC and speed up solver
PR = zeros(size(OLa.X,1), 1);

La.X = L;
% remove empty voxels touching occupied ones
Lai = fsparse(extdof,extdof,1,size(La.X)); % pressure BC is thus at oxy BC
La.X = La.X-Lai*La.X+Lai;
[La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

%% Time loop
while tt <= tspan(end)
     %% Initialise U and classify the DOFs
     adof = find(U >= cutoff_bdof); % all filled voxels
     Idof =  (U <= cutoff_bdof); % all cells with 'small' density
     
     % find tumor boundary
     nIdof = (N*Idof & ~Idof); % neighbours of Idof that are not Idof...
     Idof = (N*nIdof & Idof);  % ... neighbours of above that are Idof.
     
     idof1 = find(Idof & ~VU);         % "external" OBC1
     %idof2 = find(Idof & VU);         % "internal" OBC2
     idof = find(Idof);

     % "All DOFs" = adof + idof, like the "hull of adof"
     Adof = union(adof, idof);
     % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
     % matrix. Determine also a local enumeration, eg. [1 2 3
     % ... numel(Adof)].

     Adof_ = (1:numel(Adof))';
     [idof1_,idof_,adof_,] = ...          
        map(Adof_,Adof,idof1,idof,adof);

     %% Calculate Pressure and Oxygen systems
    
     % initial guess for the oxygen solver
     if tt == tspan(1) 
             % Calculate new oxygen distribution explicitly
             Oxy = full(fsparse([extdof; adof],1, ...
                 [ones(size(extdof)); ...
                 -lambda*full(U(adof))], ...
                 [size(OLa.X,1) 1]));
             Oxy(ndof) = 0.*Oxy(ndof);  % no consumption for necrotic cells
             Oxy(U <= cutoff_bdof) = 0; % no consumption outside
             Oxy(extdof) = 1;           % impose the BC
             Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\(M*Oxy)));
             Oxy = solve_oxygen(iOLa, M, Oxy, U, lambda, extdof, ...
                 cutoff_bdof, kappa_death, dtau, dt);
     end
     
     % Find oxygen distribution iteratively
     Oxy = solve_oxygen(iOLa, M, Oxy, U, lambda, extdof, ...
            cutoff_bdof, kappa_death, dtau, Tend/n_measurements);
     
     % find proliferating cells & dying cells (for pressure distr.)
     % non-boundary prolif. dofs
     pdof = find(U > cutoff_bdof & (Oxy > kappa_prol)); 
     ndof = find(U > cutoff_bdof & (Oxy < kappa_death)); % 'necrotic' dof
     
     % Do mapping again (could optimize this somehow?)
     [pdof_,ndof_] = ...          
     map(Adof_,Adof,pdof,ndof);
     
     % evaluate boundary curvature
     % get contour line
     [X,Y,Z] = meshdata(P(1,:),P(2,:),U');
     [cpos, ~] = contour(X, Y, Z, [cutoff_bdof,cutoff_bdof]);    
     
     % get all contour segments
     boundaries = cell(1,0);
     cn = [];
     idx = 1;                        % starting at index 1
     while idx < length(cpos)
         cn = [cn cpos(2,idx)];
         % save contour points
         boundaries{end+1} = [cpos(1,idx+1:idx+cn(end)); ...
                              cpos(2,idx+1:idx+cn(end))];
         idx = idx + cn(end) + 1;    % step to next contour line
     end
     
     % get curvature at idofs
     desired_coords = [P(1,idof); P(2,idof)]; 
     curvature = curvature2D(boundaries, desired_coords, curv_cutoff, ...
         curv_thresh, sgf_wind,sgf_deg);
     
     % get main boundary for evaluating small perturbations
     cidx = find(cn == max(cn));
     boundary_main = boundaries{cidx};  % main contour line
     
     % set Dirichlet only at idofs that are from large enough contours
     idof_keep = find(~isnan(curvature));
     curvature = curvature(idof_keep);
     
     % pressure Laplacian 
     % note: idofs are well-connected by construction...
     La.X = L(Adof,Adof);
     % ...but rescale all idof hats for Neumann condition there
     Lai2 = fsparse(idof_,idof_,1,size(La.X));  
     % Scale Laplacian on boundary
     La.X = La.X - diag(sum(Lai2*La.X,2));
     % remove empty voxels touching occupied ones, using external idof only
     Lai = fsparse(idof_(idof_keep),idof_(idof_keep),1,size(La.X)); 
     La.X = La.X-Lai*La.X+Lai;
 
     [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector'); 
     
     % Calculate new pressure distribution explicitly
     % using lumped mass matrix for the piecewise constant source-terms
     % RHS source term proportional to the over-occupancy and BCs
     Pr = full(fsparse(pdof_,1,mu_prol/D.*dM(pdof) ,...
         [size(La.X,1) 1]));     % RHS first...
     
     Pr(idof_(idof_keep)) = +sigma.*curvature;   % (set BC,
     Pr(ndof_) = -mu_death/D.*dM(ndof);          % sinks for dying dofs)
     Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\(Pr))); % ..then the solution
     
     % Extract the Adof pressures
     PR = zeros(size(OLa.X,1), 1);
     PR(Adof) = Pr;
     
     %% get move-'rates' (finite-volume approach)
     
     % Finite-volume approach - akin to the Gillespie-rates in
     % AvascularTumor.m
     [ii,jj_] = find(L); % find tumor-domain Neumann neighbours
     % remove flow that occurs outside tumor
     keep = find( (U(jj_) >= cutoff_bdof) | Idof(jj_));
     ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
     keep = find( (U(ii) >= cutoff_bdof) | Idof(ii));
     ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
        
        
     % Split approach (essentially Donor Cell Upwind method)
     % h^2 only for square dual!
     grad_loss = fsparse(ii,1, ...
         min((PR(jj_)-PR(ii)).*full(U(ii))/(hmax^2), 0), ...  % outgoing
                    Nnodes);
     grad_gain = fsparse(ii,1, ...
         max((PR(jj_)-PR(ii)).*full(U(jj_))/(hmax^2), 0), ...  %incoming
                    Nnodes);

     grad = grad_loss + grad_gain;
      
     mover = full(D*gradquotient*grad); 
     % mass growth/loss
     massr = full(fsparse(pdof,1,mu_prol.*full(U(pdof)) ,...
                                        [Nnodes 1]));
     massr(ndof) = -mu_death.*U(ndof);
     
%%  Calculate timestep dt 

     % CFL condition
     dt1 = 0.1 / max(abs(mover) + abs(massr));% ~CFL if D*gradquotient = 1 
     dt2 = hmax;    % ensures rough stability approximations for FV scheme
     dt = min(dt1,dt2);
     DT(i) = dt;    % save the time steps
     
     % end simulation if system equilibrium reached
     if dt == Inf
         warning('Equilibrium reached, dt == Inf')
         break;
     end
     
     %%  Report back and save time series of current states 
     if tspan(i+1) < tt+dt
         iend = i+find(tspan(i+1:end) < tt+dt,1,'last'); 
         
         % save relevant values 
         Usave(i+1:iend) = {U};
         
         Adofsave(i+1:iend) = {Adof};
         Prsave(i+1:iend) = {Pr};

         Oxysave(i+1:iend) = {Oxy};
         
         % get tumor boundary contours and their curvature
         % (cn, cpos, boundary, evaluated at each time-step)
         cr = sqrt(boundary_main(1,:).^2 + boundary_main(2,:).^2);
         
         % estimate contour roundness from its corresponding polygon
         pgon = polyshape(boundary_main(1,:),  boundary_main(2,:));
         perim = perimeter(pgon);
         parea = area(pgon);
         roundness(i+1:iend) = perim^2/(4*pi*parea); % = 1: perfect circle
         
         % save the Fourier transform of the contour
         Yfft = fft(cr);
         Fs = numel(cr);
         n = Fs;
         P2 = abs(Yfft/n);
         n2 = floor(n/2); % avoid fractional indexing
         P1 = P2(:, 1:n2 + 1);
         P1(:, 2:end-1) = 2*P1(:, 2:end-1);
         BFTsave(i+1:iend) = {[0:1:(n2-1);P1(1:n2)]};
         
         % simple estimate of regional sizes
         ZZ = zeros(Nnodes,1);
         ZZ(U > cutoff_bdof) = 1; % all cells part of tumor
         voltot = sum(M*ZZ);      % total volume
         ZZ(:) = 0;
         ZZ(pdof) = 1;
         volp = sum(M*ZZ);        % volume, proliferating region
         ZZ(:) = 0;
         ZZ(ndof) = 1;            % volume, necrotic region
         voln = sum(M*ZZ);
         rp(i+1:iend) = sqrt( voltot /pi);
         rq(i+1:iend) = sqrt( (voltot - volp) /pi);
         rn(i+1:iend) = sqrt( voln /pi);
         
         i = iend;
     end

     %%  Euler forward step
     
     % update densities given move rates and mass change rates
     dU = dt.*mover + dt.*massr;
     U = U + dU;
     % add artificial noise
     U = U + noise.*(sqrt(abs(dU))).*normrnd(0,1,numel(U),1);
    
     U(extdof) = 0; % avoids severe issues on comp. boundary

     %% Step in time
     
     tt = tt+dt; 
     report(tt,U,'');

     % update the visited sites
     VU = VU | U;
     
     % End sim if highly unstable or growth reaches domain edge
     if max(U) > 100 | intersect(adof, extdof)
         warning('simulation aborted due to possible instability')
         break;
     end
end

%report(tt,U,'done');

%return;

%% Plots

if ~exist('cutoff_bdof', 'var')
    % load parameters
    lambda = all_parameters(1);
    kappa_prol = all_parameters(2);
    mu_prol = all_parameters(3);
    kappa_death = all_parameters(4);
    mu_death = all_parameters(5);
    cutoff_bdof = all_parameters(6);
    sigma = all_parameters(7);
end

if ~exist('tstamps', 'var')
    % set default times for figures
    tstamps = [2 round(0.8*numel(tspan)) numel(tspan)];
end

i = 1; % figure index
% plot solution at times tstamps
for iii = tstamps
    fig=figure(i);
    set(gcf,'Position',[100 100 340/3 240/2]);
    clf,

    Umat=full(cell2mat(Usave));

    % background
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
     'EdgeColor','none');

    hold on,
    axis(0.75.*[-1 1 -1 1]); axis square, axis off

    % colour tumor volumes after growth/death/quiescence

    pdof = find(Usave{iii} > cutoff_bdof & ... % non-boundary prolif. dofs
     (Oxysave{iii} > kappa_prol));  
    qdof = find(Usave{iii} > cutoff_bdof & ... % quiescent dofs
     (Oxysave{iii} <= kappa_prol) & (Oxysave{iii} >= kappa_death));
    ndof = find(Usave{iii} > cutoff_bdof & ... % 'necrotic' dof
     (Oxysave{iii} < kappa_death)); 

    % nonzero density outside tumor
    outdof = find(Usave{iii} & Usave{iii} <= cutoff_bdof); 

    % plot regions
    c = graphics_color('vermillion');
    patch('Faces',R(pdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');     
    c = graphics_color('bluish green');
    patch('Faces',R(qdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0 0 0];
    patch('Faces',R(ndof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0.6 0.6 0.6];   % darker grey than background
    patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 

    % mark domain origo by crosshairs
    plot([0 0], [-1 1],'w')
    plot([-1 1], [0 0] ,'w')
    drawnow;
    
    set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
    i = i+1;
end
 
% plot volumetric growth
PDE = true;
RegionalCharacteristicsVisualiser
 
return

%% Save data & animate
saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'Prsave', ...
                  {Prsave}, 'Oxysave', {Oxysave}, 'tspan', {tspan}, ...
                  'R', {R}, 'V', {V}, 'N', {N}, 'Pr', {Pr}, 'Adof', ...
                  {Adof},'adof', {adof}, 'adof_', {adof_}, 'idof', ...
                  {idof},'idof_', {idof_}, 'Nnodes',{Nnodes}, 'P', ...
                  {P}, 'noise', {noise}, 'rng_seed', {seed}, 'rp', ...
                  {rp}, 'rq', {rq}, 'rn', {rn}, 'k', {k}, 'hmax', ...
                  {hmax}, 'Tend', {Tend}, 'all_parameters', {allpar}, ...
                  'roundness', {roundness}, 'radius', {radius}, ...
                  'BFTsave', {BFTsave}, 'tstamps', {tstamps});
filename = "PDEtumor_test";
filename_saveData = filename + ".mat";
save(filename_saveData,'-struct','saveData');

return;
% animate tumor growth
animate_growth(filename_saveData, "PDE", false)
