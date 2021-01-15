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
%   Permeability: Drate1 describes the rate diffusion rate of tumour
%   cells invading previously unvisited voxels. Drate2 is the rate
%   cells move into previously occupied but currently empty
%   voxels. Drate3 is the rate cells move into voxels that are already
%   occupied.

% M.C. Jayaweera & A. Graf Brolund 2020-12(revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

clear;
clc;
close all;

% profile on

%% Initial experiment setup
% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
mesh_type = 1;
[P,E,T,gradquotient] = basic_mesh2(mesh_type,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

D = 1; % D rate, the rate of which cells move in the domain 
Drate_ = [0.01;25];

% simulation interval
Tend = 100;
tspan = linspace(0,Tend,101);
timescaling=0.005;
report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

% boundary conditions
alpha = 1e+4; % weighting parameter for Robin BC
alpha_inv = 1/alpha; % inverse is the value used in simulation

% Set which experiment to use
exp = 0;

%Experiments
if exp == 0
    % Normal run------------------------------------
    start_value = 1;
    radius = 0.05;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
    
    % Initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]);      %Intialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]);   %Initialize
    
elseif exp ==1
    % Initial: relaxation---------------------------
    start_value = 10;
    radius = 0.07;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0;        % rate of death
    r_degrade = 0;     % rate of degradation for already dead cells
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
elseif exp==2 
    %initial: simulation---------------------------
    start_value = 1;
    radius = 0.15;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells 
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    
elseif exp==3 
    %initial: Quick, concentrated------------------
    start_value = 1;
    radius = 0.09;
    
    %Tend=30;
    cons = 0.15;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.3;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.3;        % rate of death
    r_degrade = 0.05;     % rate of degradation for already dead cells 
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
end

cutoff_bdof = 0.1;
cutoff_deg = 0.0001;
cutoff_remain = 0.01;

%% Other simulation setup
% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);       %N gives the neighbours
Mgamma = assemble_Mgamma(P,T);      % Assebmle boundary matrix
Mgamma = Mgamma./dM;
neigh_Mgamma = (Mgamma ~= 0) - speye(size(Mgamma));
neigh = full(sum(N,2));
N_vec = zeros(size(N,1),4);         % Matrix used to find sdof_m
for k = 1:size(N,1)
    temp = find(N(k,:));
    if length(temp) == 2
        temp = [temp, temp(1),temp(2)];
    elseif length(temp) == 3
        temp = [temp, temp(1)];
    end
    N_vec(k,:) = temp;
end

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;
Udsave = cell(1,numel(tspan));
Udsave{1} = U_dead;

Oxysave = cell(1,numel(tspan));
bdofsave = cell(1,numel(tspan));
sdofsave = cell(1,numel(tspan));
sdofbsave = cell(1,numel(tspan));

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);
inspect_rate_toIdof = cell(3,numel(tspan));
% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   

[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

%% Time loop
while tt <= tspan(end)
    %% Init U and U_dead and classify the DOFs
    U = U_new;
    U_dead = U_deadnew;
    U_and_U_dead = U | U_dead;
    
    adof = find(U_and_U_dead); % all filled voxels 
    % singularly occupied voxels on the boundary: 
    bdof_m = find(N*(U_and_U_dead ~= 0) < neigh & (U > cutoff_bdof & ...
        U <= 1));
    sdof = find(U > 1); %  sdof on the boundary
    sdof_b = find(N*(U_and_U_dead ~=0) < neigh & (U > 1));
    % voxels with more than concentration 1 in them which may move, 
    % with a voxel containing less number of cells next to it:
%     sdof_m = find(sum(N.*(U')<(U)&(N&N),2).*(U>1));
    sdof_m = find(U - min(U(N_vec),[],2) > 0 & U>1);
    Idof = (N*(U_and_U_dead ~= 0) > 0 & U_and_U_dead == 0); % empty voxels touching occupied ones
    idof1 = find(Idof & ~VU);         % "external" OBC1
    idof2 = find(Idof & VU);          % "internal" OBC2
    idof3 = find(~VU & neigh_Mgamma*VU > 0); % boundary around "visited voxels"
    idof3 = setdiff(idof3,idof1);
    idof = find(Idof);
    ddof = find(U_dead > 0);   %degrading voxels
    
    % Add idof3
    idof1 = [idof1;idof3];
    idof = [idof;idof3];
    
    % "All DOFs" = adof + idof, like the "hull of adof"
    Adof = [adof; idof];
    
    % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
    % matrix. Determine also a local enumeration, eg. [1 2 3
    % ... numel(Adof)].
    Adof_ = (1:numel(Adof))';
    [bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_, sdof_b_,ddof_] = ...
        map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof,adof,sdof_b,ddof);
    
    %% Calculate Pressure and Oxygen systems
    %Pressure and oxygen calculation
%     if ~isempty(idof2) 
%         disp("JAH!");
%         pause
%     end
    if updLU
        % pressure Laplacian
        La.X = L(Adof,Adof);
        %%% ADD BC to LHS
        % Robin to external idof (idof1)
        Lai1 = fsparse(idof1_,idof1_,1,size(La.X));
%         Lai2 = fsparse(idof2_,idof2_,1,size(La.X));
        a_Lai1 = speye(size(Lai1)) - Lai1;

        % Scale Laplacian on boundary
        La.X = La.X - Lai1.*La.X;
        La.X = La.X - diag(sum(Lai1*La.X,2));
%         La.X = La.X - Lai2*La.X + Lai2;

        % Get local Mgamma for all active dofs
        Mgamma_b = Mgamma(Adof,Adof);

        % Get only idof1 part of Mgamma_b and set the diagonal as the sum of
        % all non-diagonal elements times 2
        Mgamma_b = Lai1*Mgamma_b*Lai1 - Mgamma_b.*Lai1;
        Mgamma_b = Mgamma_b + diag(2*sum(Mgamma_b,2));

        % Put together the LHS and remove the connection from adof to idof1
        % (this sets Dirichlet for adof but keeps the Robin for the boundary)
        La.X = La.X + alpha_inv*Mgamma_b - a_Lai1*La.X*Lai1;
        % LU-factorize
        [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
    end

    % RHS source term proportional to the over-occupancy and BCs
    Pr = full(fsparse(sdof_,1,(U(sdof)-1)./dM(sdof), ...     % Take U_dead into consideration?
        [size(La.X,1) 1]));     % RHS first...
    Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

    % RHS source term proportional to the over-occupancy and BCs
    Oxy = full(fsparse([extdof; adof],1, ...
        [ones(size(extdof)); ...
        -cons*full(U(adof)./dM(adof))], ... 
        [size(OLa.X,1) 1]));
    Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));

    if isempty(Oxysave{1})
        Oxysave{1}=Oxy;
    end
    %% Move calculations
    %sdof_m  
    rates_sdof = zeros(length(Adof),1);
    [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
    Pr_diff__ = max(Pr(sdof_m_(ii))-Pr(jj_),0);    %proportional to over-occupancy
    grad_sdof = fsparse(ii,1,Pr_diff__.*Drate_(1+VU(Adof(jj_))), numel(sdof_m));
    rates_sdof(sdof_m_) = -gradquotient*grad_sdof;
    grad_N = fsparse(jj_,1, Pr_diff__.*Drate_(1+VU(Adof(jj_))), numel(Adof));
    rates_sdof = rates_sdof + gradquotient*grad_N;
    
    % bdof_m
    rates_bdof = zeros(length(Adof),1);
    %check if boundary is updated. If no bdofs exist, skip bdof calculation
    if (sum(bdof_m) + sum(sdof_b)) ~= 0
        for ind=1:length(bdof_m_)
            % movement of a boundary (singly occupied) voxel
            ix = bdof_m(ind);
            ix_ = bdof_m_(ind);

            jx_ = find(N(ix,Adof));
            % (will only move into an empty voxel:)
            jx_ = jx_(U_and_U_dead(Adof(jx_)) == 0);
            Pr_diff = max(Pr(ix_)-Pr(jx_),0).*Drate_(1+VU(Adof(jx_)));%*U(ix)*;    %proportionellt mot over-occupancy

            rates_bdof(jx_) = rates_bdof(jx_) + Pr_diff;
            rates_bdof(ix_) = -sum(Pr_diff); 
        end
    end
    
%     if sum(rates_bdof(idof_)) > 0
%         disp("Hold up");
%     end
    
    %% Change calculation
    %proliferation-----------------------------
    ind_prol = find((Oxy > cutoff_prol));
    prol_conc = r_prol*U(ind_prol);
    
    
    %death--------------------------------------
    ind_die = find(Oxy < cutoff_die); %index for dying cells
    dead_conc = r_die*U(ind_die);
        
    %degradation--------------------------------
    degrade_conc = U_deadnew(ddof)*r_degrade; 
    
    %% Intensity calculation
    % bdof_m
    moveb = rates_bdof(bdof_m_);
    % sdof_m
    moves = rates_sdof(sdof_m_);
    % birth
    birth = (ind_prol>0)*r_prol; 
    % death
    death = (ind_die>0)*r_die; 
    % degradation
    degrade = (r_degrade*(U_deadnew(ddof)>0));
    % total
    intens = [birth; death; degrade; moves; moveb];
    %%%#########
    intens_print__ = [sum(birth); sum(death); sum(degrade); sum(moves); sum(moveb)]; 
    %% Calculate timestep dt
    ind_rates_sdof_n = find(rates_sdof(sdof_m_)<0);
    ind_rates_bdof_n = find(rates_bdof(bdof_m_)<0);

    dt_death = U_new(ind_die)./(dead_conc);
    dt_sdof = U_new(sdof_m(ind_rates_sdof_n))./(-rates_sdof(sdof_m_(ind_rates_sdof_n)));
    dt_bdof = U_new(bdof_m(ind_rates_bdof_n))./(-rates_bdof(bdof_m_(ind_rates_bdof_n)));
    
    if (sum(isinf(dt_sdof))>0)
        inf;
    end
    dt_test = (min([dt_death; dt_sdof; dt_bdof;(0.1*Tend)]));
    dt = dt_test*timescaling;
    %% Save time series of current step
    % Report back and save time series 
    if tspan(i+1) < tt+dt
        iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
        
        % save relevant values 
        Usave(i+1:iend) = {U};
        Udsave(i+1:iend) = {U_dead};

        Oxysave(i+1:iend) = {Oxy};

        bdofsave(i+1:iend) = {bdof_m};
        sdofsave(i+1:iend) = {sdof_m};
        sdofbsave{i+1:iend} = {sdof_b};
        
        % the rates
        inspect_rate_toIdof(:,i) = {rates_bdof;sdof_m_;idof_};
        
        i = iend;
        
        % monitor the maximum outlier cell:
        max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));
    end
    
    %% Make the Euler forward step
    %Proliferation
    U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;

    %Death
    U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
    U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
    ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
    U_new(ind_cutoff) = 0;

    % Degradation
    U_deadnew(ddof) = U_deadnew(ddof) - degrade_conc*dt;
    U_deadnew(U_deadnew < cutoff_deg) = 0; % remove cells below cutoff_deg

    % sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt;

    % bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt;
    %% Step in time
    tt = tt+dt;
    report(tt,U,'');
    
    % update the visited sites
    VU = VU | U_new;
end
report(tt,U,'done');

%%
TumorGraphicsResult;

%% SAVE DATA
saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'tspan', {tspan}, ...
    'R', {R}, 'V', {V}, 'N', {N}, 'Udsave', {Udsave}, ...
    'max_radius', {max_radius}, ...
    'alpha', {alpha}, 'Pr', {Pr}, 'Adof', {Adof},  ...
    'adof', {adof}, 'adof_', {adof_}, 'idof', {idof}, 'idof_', {idof_}, ...
    'Nvoxels',{Nvoxels},'inspect_rates',{inspect_rate_toIdof}, ...
    'P', {P}, 'bdof_m', {bdof_m}, 'bdof_m_', {bdof_m_}, ...
    'sdof_m', {sdof_m},'sdof_m_', {sdof_m_}, 'gradquotient', {gradquotient}, ...
    'Tend', {Tend});
filename_ = "alpha" + erase(sprintf('%0.0e',alpha),'.');
filename_ = filename_ + "_" + strjoin(string(fix(clock)),'-');
filename_saveData = "saveData_" + filename_ + ".mat";
save(filename_saveData,'-struct','saveData');

return;