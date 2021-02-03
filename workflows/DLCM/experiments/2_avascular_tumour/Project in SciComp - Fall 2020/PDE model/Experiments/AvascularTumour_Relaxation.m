% Simulation of an avascular tumour model.
%
%   Avascular tumour growth, relaxation experiment: An initial circular 
%   population of cells lie crowded together (concentration>1). When 
%   released, the pressure exerted between the cells will have them move 
%   out in the domain.
%
%   No proliferation, death, or degradation is allowed.  

% C. Jayaweera & A. Graf Brolund 2021-01 (revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

clear;
clc;
close all;

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for... 
[V,R] = mesh2dual(P,E,T,'voronoi');             %cartesian mesh

D = 1; % D_rate, the rate with which cells move in the domain. 
        % currently the rate is the same for visited voxels and non-visited

% simulation interval
Tend = 100;
tspan = linspace(0,Tend,101);
timescaling = 0.005;

% initial population: circular blob of living cells    
start_value = 10; % cell concentrations in the initial blob
radius = 0.07;
r = sqrt(P(1,:).^2+P(2,:).^2);
ii = find(r < radius); % radius of the initial blob
U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]); % initialize
U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); 
U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); 

% parameters
cutoff_bdof = 0.1;    % lower bound for bdof

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);       %N gives the neighbours 
neigh = full(sum(N,2));

N_vec = zeros(size(N,1),4);          % neighbour matrix used to find sdof_m
for k = 1:size(N,1)
    temp = find(N(k,:));
    if length(temp) == 2
        temp = [temp, temp(1),temp(2)];
    elseif length(temp) == 3
        temp = [temp, temp(1)];
    end
    N_vec(k,:) = temp;
end

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;
Udsave = cell(1,numel(tspan));
Udsave{1} = U_dead;

% keeps track of dofs for figures and understanding
bdofsave = cell(1,numel(tspan));
sdofsave = cell(1,numel(tspan));
sdofbsave = cell(1,numel(tspan));

tt = tspan(1);
i = 1;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);

while tt <= tspan(end)
    U = U_new;
    U_dead = U_deadnew;
    
    %% Init U and U_dead and classify the DOFs
    U = U_new;
    U_dead = U_deadnew;
    U_and_U_dead = U | U_dead;
  
    %Classification of the DOFs
    
    adof = find(U_and_U_dead); % all filled voxels 
    % singularly occupied voxels on the boundary: 
    bdof_m = find(N*(U_and_U_dead ~= 0) < neigh & (U > cutoff_bdof & ...
        U <= 1));
    sdof = find(U > 1); % source voxels,concentration more than 1
    % sdof on the boundary
    sdof_b = find(N*(U_and_U_dead ~=0) < neigh & (U > 1)); 
    % voxels with more than concentration 1 in them which may move, 
    % with a voxel containing a lower concentrations next to it:
    sdof_m = find(U - min(U(N_vec),[],2) > 0 & U>1);
    % empty voxels touching occupied ones  
    Idof = (N*(U_and_U_dead ~= 0) > 0 & U_and_U_dead == 0);      
    idof1 = find(Idof & ~VU);         % "external" OBC1
    idof2 = find(Idof & VU);          % "internal" OBC2
    idof = find(Idof);
    
    % "All DOFs" = adof + idof, like the "hull of adof"
    Adof = [adof; idof];
    % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
    % matrix. Determine also a local enumeration, eg. [1 2 3
    % ... numel(Adof)].

    Adof_ = (1:numel(Adof))';
    [bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_, sdof_b_] = ...          
      map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof,adof,sdof_b);
    
    %% Calculate Pressure

    % pressure Laplacian
    La.X = L(Adof,Adof);
    %remove emtpy voxels touching occupied ones
    Lai = fsparse(idof_,idof_,1,size(La.X)); 
    La.X = La.X-Lai*La.X+Lai;
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

    % RHS source term proportional to the over-occupancy and BCs
    Pr = full(fsparse(sdof_,1,(U(sdof)-1)./dM(sdof),... %equilibrium at U=1
        [size(La.X,1) 1]));   % RHS first...                           
    Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

    %% Movement calculations
    
    % movement of cells in sources, sdof_m
    rates_sdof = zeros(length(Adof),1);
    [ii,jj_] = find(N(sdof_m,Adof)); % neighbours
    % pressure difference between cells and its neighbours    
    Pr_diff__ = max(Pr(sdof_m_(ii))-Pr(jj_),0);   
    grad_sdof = fsparse(ii,1,Pr_diff__*D, numel(sdof_m));
    rates_sdof(sdof_m_) = -gradquotient*grad_sdof; % sources lose cells
    grad_N = fsparse(jj_,1, Pr_diff__*D, numel(Adof));
    rates_sdof = rates_sdof + gradquotient*grad_N; % neighbours gain cells
                         
    % movement of cells on the boundary, not over-occupied, bdof_m
    % this loop is similar to what is done above for sdof and can be 
    % rewritten in a similar fashion. Since this calculation is less 
    % common, it is not a big concern for efficiency
    rates_bdof = zeros(length(Adof),1);
    for ind=1:length(bdof_m_)
        ix = bdof_m(ind);
        ix_ = bdof_m_(ind);

        jx_ = find(N(ix,Adof)); % neighbour to the voxel
        jx_ = jx_(U_and_U_dead(Adof(jx_)) == 0); % empty neighbours
        Pr_diff = max(Pr(ix_)-Pr(jx_),0); % pressure difference

        rates_bdof(ix_) = -sum(D*Pr_diff); % loses cells        
        rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff; % gain cells
    end
 
    %%  Calculate timestep dt 
    
    ind_rates_sdof_n = find(rates_sdof(sdof_m_)<0); % affected voxels
    ind_rates_bdof_n = find(rates_bdof(bdof_m_)<0);

    % find the largest possible time step while avoiding U<0
    dt_sdof = U_new(sdof_m(ind_rates_sdof_n))./ ... 
        (-rates_sdof(sdof_m_(ind_rates_sdof_n)));
    dt_bdof = U_new(bdof_m(ind_rates_bdof_n))./ ... 
        (-rates_bdof(bdof_m_(ind_rates_bdof_n)));

    dt = min([dt_sdof; dt_bdof;(0.1*Tend)])*timescaling; % scale dt smaller
   
    %% Report back and save time series of current states 

    if tspan(i+1) < tt+dt
        iend = i+find(tspan(i+1:end) < tt+dt,1,'last');

        % save relevant values 
        Usave(i+1:iend) = {U};
        Udsave(i+1:iend) = {U_dead};

        bdofsave(i+1:iend) = {bdof_m};
        sdofsave(i+1:iend) = {sdof_m};
        sdofbsave{i+1:iend} = {sdof_b};

        i = iend;
    end
    
    %% Euler forward step
    
    % movement of cells in sources, sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt;

    % movement of cells in boundary voxels, bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt;
    
    tt = tt+dt;
%     report(tt,U,'');
    
    % update the visited sites
    VU = VU | U;
end
% report(tt,U,'done');

%% Create a GIF animation
Mnormal = struct('cdata',{},'colormap',{});
figure(1), clf,

Umat=full(cell2mat(Usave));
colorbar
caxis([0 max(max(Umat))])
colorlabel('Concentration of cells, U')
for i = 1:numel(Usave)
    % background
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    % colour living voxels after concentration level
    ii = find(Usave{i}>0);
    c = Umat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c, ... 
        'FaceColor','flat');     
    
    % color (fully) dead voxels black
    ii = find(Usave{i} == 0 & Udsave{i} > 0);
    p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);
    legend(p_dead,'dead')
    
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
    Mnormal(i) = getframe(gcf);
end

% save the GIF
movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'Tumour.gif', ...
          'delaytime',0.1,'loopcount',0);
