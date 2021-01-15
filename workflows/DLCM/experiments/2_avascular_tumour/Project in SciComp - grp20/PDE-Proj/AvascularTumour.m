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

% M.C. Jayaweera & A. Graf Brolund 2021-01(revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

clear;
clc;
close all;

%profile on       %Profiler to check which functions are most
                  %computationally heavy

%% Initial experiment setup
% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

D=1; %D_rate, the rate of which cells move in the domain  

% Simulation interval
Tend = 11;                  %Final time step
tspan = linspace(0,Tend,101);
timescaling=0.005;          %Time scaling

report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

%     The user specified cutoff and rate parameters for the proliferation,
%     death, degradation and consumption rules.
start_value = 1;
radius = 0.05;

cons = 0.0015;        % consumption of oxygen by cells
cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
r_prol = 0.125;       % rate of proliferation of singly occupied voxels
cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
r_die = 0.125;        % rate of death
r_degrade = 0.01;     % rate of degradation for already dead cells

% cutoff parameters
cutoff_bdof = 0.1;
cutoff_deg = 0.0001;
cutoff_remain = 0.01;

% Initial population: circular blob of living cells
r = sqrt(P(1,:).^2+P(2,:).^2);
ii = find(r < radius); % radius of the initial blob
U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]);      %Intialize
U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]);   %Initialize

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary
    
% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);       %N gives the neighbours 
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

%%
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

%% Time loop
while tt <= tspan(end)
    %% Init U and U_dead and classify the DOFs
    U = U_new;
    U_dead = U_deadnew;
    U_and_U_dead = U | U_dead;
  
    %Classification of the DOFs
    adof = find(U_and_U_dead); % all filled voxels 
    % singularly occupied voxels on the boundary: 
    bdof_m = find(N*(U_and_U_dead ~= 0) < neigh & (U > cutoff_bdof & ...
        U <= 1));
    sdof = find(U > 1); %  sdof on the boundary
    sdof_b = find(N*(U_and_U_dead ~=0) < neigh & (U > 1));
    % voxels with more than concentration 1 in them which may move, 
    % with a voxel containing less number of cells next to it:
    sdof_m = find(U - min(U(N_vec),[],2) > 0 & U>1);
    Idof = (N*(U_and_U_dead ~= 0) > 0 & U_and_U_dead == 0); % empty voxels touching occupied ones idof1 = find(Idof & ~VU);         % "external" OBC1
    idof1 = find(Idof & ~VU);         % "external" OBC1
    idof2 = find(Idof & VU);          % "internal" OBC2
    idof = find(Idof);
    ddof = find(U_dead > 0);            %degrading voxels
    
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
    
    % pressure Laplacian
    La.X = L(Adof,Adof);
    Lai = fsparse(idof_,idof_,1,size(La.X)); %remove emtpy voxels touching occupied ones
    La.X = La.X-Lai*La.X+Lai;
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

    % RHS source term proportional to the over-occupancy and BCs
    Pr = full(fsparse(sdof_,1,(U(sdof)-1)./dM(sdof), ...             %the equilibrium value of pressure is U=1 or U(sdof-1)
        [size(La.X,1) 1]));     % RHS first...
    Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

    % RHS source term proportional to the over-occupancy and BCs
    Oxy = full(fsparse([extdof; adof],1, ...
        [ones(size(extdof)); ...
        -cons*full(U(adof)./dM(adof))], ... 
        [size(OLa.X,1) 1]));
    Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));

%% Move calculations
    %sdof and bdof movement calculations
    
    %sdof_m  
    rates_sdof = zeros(length(Adof),1);
    [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
    Pr_diff__ = max(Pr(sdof_m_(ii))-Pr(jj_),0);    %proportional to over-occupancy
    grad_sdof = fsparse(ii,1,Pr_diff__*D, numel(sdof_m));
    rates_sdof(sdof_m_) = -gradquotient*grad_sdof;
    grad_N = fsparse(jj_,1, Pr_diff__*D, numel(Adof));
    rates_sdof = rates_sdof + gradquotient*grad_N;
    
%     % bdof_m
%     rates_bdof = zeros(length(Adof),1);
% 
%     %check if boundary is updated. If no bdofs exist, skip bdof calculation
%     if sum(bdof_m) == 0
%         updLU = true;
%         return          %!Check this expression if integrated into larger code
%     else
%        
    % bdof_m
    rates_bdof = zeros(length(Adof),1);
    %check if boundary is updated. If no bdofs exist, skip bdof calculation
%         for ind=1:length(bdof_m_)
%             % movement of a boundary (singly occupied) voxel
%             ix = bdof_m(ind);
%             ix_ = bdof_m_(ind);
% 
%             jx_ = find(N(ix,Adof));
%             % (will only move into an empty voxel:)
%             jx_ = jx_(U_and_U_dead(Adof(jx_)) == 0);
%             Pr_diff = max(Pr(ix_)-Pr(jx_),0);
% 
%             rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff;
%             rates_bdof(ix_) = -sum(D*Pr_diff); 
%         end
        
        for ind=1:length(bdof_m_)
            % movement of a boundary (singly occupied) voxel
            ix = bdof_m(ind);
            ix_ = bdof_m_(ind);

            jx_ = find(N(ix,Adof));
            % (will only move into an empty voxel:)
            jx_ = jx_(U_and_U_dead(Adof(jx_)) == 0);
            %jx_ = jx_(U(Adof(jx_)) == 0 & U_dead(Adof(jx_)) == 0);
            Pr_diff = max(Pr(ix_)-Pr(jx_),0);%*U(ix);    %proportionellt mot over-occupancy

            rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff;
            rates_bdof(ix_) = -sum(D*Pr_diff); 
        end

    
        
%         rates_bdof = zeros(length(Adof),1);
%         [ii,jj_] = find(N(bdof_m,Adof)); % neighbours...
%         %keep = find(U(Adof(jj_)) < 1);   % ...to move to
%         %ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
%         % remove any possibly remaining negative rates
%         %Pr_diff__ = max(Pr(bdof_m_(ii))-Pr(jj_),0).*(U(bdof_m(ii))-1);    %proportionellt mot over-occupancy
%         Pr_diff__ = max(Pr(bdof_m_(ii))-Pr(jj_),0);    %proportionellt mot over-occupancy
% 
%         grad_bdof = fsparse(ii,1,Pr_diff__*D, numel(bdof_m)); 
%         %moves = full(gradquotient*grad);
%         rates_bdof(bdof_m_) = -gradquotient*grad_bdof;
%     
%         grad_N = fsparse(jj_,1, Pr_diff__*D, numel(Adof)); 
%         rates_bdof = rates_bdof + gradquotient*grad_N;
% %         updLU = true;
        
%     end
    updLU = true;

%%  Change calculation
    %Change calculation of proliferation, death and degradation

    %proliferation-----------------------------
    ind_prol = find((Oxy > cutoff_prol));       %index of proliferating cells
    prol_conc = r_prol*U(ind_prol);

    %death--------------------------------------
    ind_die = find(Oxy < cutoff_die);           %index for dying cells
    dead_conc = r_die*U(ind_die);

    %degradation--------------------------------
    degrade_conc = U_deadnew(ddof)*r_degrade; 

%%  Calculate timestep dt  
    % Time stsep calculation
    ind_rates_sdof_n = find(rates_sdof(sdof_m_)<0);
    ind_rates_bdof_n = find(rates_bdof(bdof_m_)<0);

    dt_death = U_new(ind_die)./(dead_conc);
    dt_sdof = U_new(sdof_m(ind_rates_sdof_n))./(-rates_sdof(sdof_m_(ind_rates_sdof_n)));
    dt_bdof = U_new(bdof_m(ind_rates_bdof_n))./(-rates_bdof(bdof_m_(ind_rates_bdof_n)));
    
    if (sum(isinf(dt_sdof))>0)
        inf;
    end
    dt_unscaled = (min([dt_death; dt_sdof; dt_bdof;(0.1*Tend)]));
    dt = dt_unscaled*timescaling;

%%  Save time series of current step
    % Report back and save time series of current states 
    if tspan(i+1) < tt+dt
        iend = i+find(tspan(i+1:end) < tt+dt,1,'last');

        % save relevant values 
        Usave(i+1:iend) = {U};
        Udsave(i+1:iend) = {U_dead};

        Oxysave(i+1:iend) = {Oxy};

        bdofsave(i+1:iend) = {bdof_m};
        sdofsave(i+1:iend) = {sdof_m};
        sdofbsave{i+1:iend} = {sdof_b};

        i = iend;

        % monitor the maximum outlier cell:
        max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));
    end

%%  Euler forward step
    %Proliferation
    U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;

    %Death
    U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
    U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
    ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));    %remove cells below cutoff_remain and cutoff_die
    U_new(ind_cutoff) = 0;

    % Degradation
    U_deadnew(ddof) = U_deadnew(ddof) - degrade_conc*dt;
    U_deadnew(U_deadnew < cutoff_deg) = 0;          %remove cells below cutoff_deg

    % sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt;

    % bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt;


%% Step in time
    tt = tt+dt;                %Increase time with time step
    report(tt,U,'');
    
    % update the visited sites
    VU = VU | U_new;
end
%report(tt,U,'done');

% return; 

%profile off                  %End of profiler
%profile report

%% Graphics

tumour_figure=1;

if tumour_figure==1
    % create a GIF animation
    Mnormal = struct('cdata',{},'colormap',{});
    fig = figure(11), 
    clf,
    Umat=full(cell2mat(Usave));
    %colorbar('southoutside')
    colorbar;
    caxis([0 max(max(Umat))])
    colorlabel('Concentration of cells, U')
    
    title=0;        %Time in gif
    snapshot =1;    %Save 5 snapshots 

    for i = 1:numel(Usave)
        %clf
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
        %patch('Faces',R,'Vertices',V,'FaceColor','none', ...
        %%transparent background
        %    'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off

        ii = find(Usave{i}>0);
        c = Umat(ii,i);
        patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');     

        ii = find(Usave{i} == 0 & Udsave{i} > 0);
        p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
        %legend(p_dead,'dead')

        if title==1
            %title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
            title(sprintf('Time = %d',tspan(i)));
        end
        drawnow;
        Mnormal(i) = getframe(gcf);
        if snapshot ==1
            %Save 5 snapshots of the tumor progression
            if i~= [1 ceil([0.24 0.49 0.74 1]*numel(Usave))]
            elseif i==1
              filename = 'T=1.eps';
              print(fig,filename,'-painters','-depsc'); 
            else
              ii=ceil(i*(Tend/numel(Usave)));
              %filename = ['T=' num2str(ii) '.eps'];
              %print(fig,filename,'-painters','-depsc'); 
               filename = ['T=' num2str(ii) '.png'];
               print(fig,filename,'-painters','-dpng'); 
            end
        end
        % saves the GIF
        %movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'Tumour.gif', ...
        %          'delaytime',0.1,'loopcount',0);
              
    end
end
%% SAVE DATA
saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'tspan', {tspan}, ...
    'R', {R}, 'V', {V}, 'N', {N}, 'Udsave', {Udsave}, ...
    'max_radius', {max_radius}, 'Pr', {Pr}, 'Adof', {Adof},  ...
    'adof', {adof}, 'adof_', {adof_}, 'idof', {idof}, 'idof_', {idof_}, ...
    'Nvoxels',{Nvoxels}, ...
    'P', {P}, 'bdof_m', {bdof_m}, 'bdof_m_', {bdof_m_}, ...
    'sdof_m', {sdof_m},'sdof_m_', {sdof_m_}, 'gradquotient', {gradquotient}, ...
    'Tend', {Tend});
%filename_ = "alpha" + erase(sprintf('%0.0e',alpha),'.');
filename_ = filename_ + "_" + strjoin(string(fix(clock)),'-');
filename_saveData = "testShit/saveData_" + filename_ + ".mat";
save(filename_saveData,'-struct','saveData');

return;
