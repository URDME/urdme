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

% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

clear;
clc;
close all;

deadfigure=0;
normalfigure=0;

% simulation interval
Tend = 50;
tspan = linspace(0,Tend,101);
% report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

%     The user specified cutoff and rate parameters for the proliferation,
%     death, degradation and consumption rules.

cutoff_bdof = 0.1;
cutoff = 0.0005;
cutoff_deg = 0.01;
cutoff_remain = 0.001;

rates_type = 1; %Relaxation
%rates_type = 2; %

if rates_type == 1
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
elseif rates_type == 2
    cons = 0;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.3;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
    
elseif rates_type == 3
    cons = 0.1;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.3;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells 
end

% Permeability parameters.
Drate1 = 0.01;     % into free matrix
Drate2 = 25;       % into already visited matrix
Drate3 = 0.01;     % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T); %N=raden ger alla voxlar, 1or för de som är grannar med radens voxel, tom eller ej
neigh = full(sum(N,2));

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

init = 1;
if init == 1
    start_value = 1;
    radie = 0.11;
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radie); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera

else
    start_value = 1;
    radie = 0.05;
    % initial population: circular blob of dead cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radie); % radius of the initial blob
    U = fsparse(ii(:),1,0,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,0,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,1,[Nvoxels^2 1]);
    U_deadnew = fsparse(ii(:),1,1,[Nvoxels^2 1]); %initiera
end

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

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
% Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   %add Dirichlet(?) oxygen at the oxygen source/outer circle
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

while tt <= tspan(end)
    U = U_new;
    U_dead = U_deadnew;
    %   Ncells=full(sum(Usave{i}~=0))
    %   Ndead=full(sum(Usave{i}==-1))
    %   Nalive=full(sum(Usave{i}>0))
    
    % classify the DOFs
    adof = find(U|U_dead); % all filled voxels (U_dead not necessary if U=-1)
    % singularly occupied voxels on the boundary: 
    bdof_m = find(N*(U ~= 0 | U_dead ~= 0) < neigh & (U > cutoff_bdof & U <= 1));%(U > 0 & U <= 1));
    sdof_b = find(N*(U ~= 0 | U_dead ~= 0) < neigh & (U > 1));
    sdof_m = intersect(find(sum(N.*U'<U & boolean(N),2)), find(U > 1)); 
    
    sdof = find(U > 1); % voxels with 2 cells
    % voxels with 2 cells in them _which may move_, with a voxel
    % containing less number of cells next to it (actually 1 or 0):
    sdof_m_old = find(N*(U > 1 | U_dead >1)<neigh & U > 1); %kanske större än 1
    Idof = (N*(U ~= 0) > 0 & U == 0); % empty voxels touching occupied ones
    idof1 = find(Idof & ~VU); % "external" OBC1
    idof2 = find(Idof & VU);  % "internal" OBC2
    idof = find(Idof);
    ddof = find(U_dead > 0);   %degrading voxels
    
    % "All DOFs" = adof + idof, like the "hull of adof"
    Adof = [adof; idof];
    % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
    % matrix. Determine also a local enumeration, eg. [1 2 3
    % ... numel(Adof)].
    Adof_ = (1:numel(Adof))';
    [bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_,sdof_b_,ddof_] = ...
        map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof,adof,sdof_b,ddof);
    
    if updLU
        % pressure Laplacian
        La.X = L(Adof,Adof);
        Lai = fsparse(idof_,idof_,1,size(La.X)); %remove emtpy voxels touching occupied ones
        La.X = La.X-Lai*La.X+Lai;
        [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
        
        updLU = false; % assume we can reuse
    end
    
    % RHS source term proportional to the over-occupancy and BCs
    Pr = full(fsparse(sdof_,1,(U(sdof)-1)./dM(sdof), ...
        [size(La.X,1) 1]));     % RHS first...
%     Pr = full(fsparse(sdof_,1,1./dM(sdof), ...
%         [size(La.X,1) 1]));     % RHS first...
    Pr_RHS = Pr;
    Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution
    
    % RHS source term proportional to the over-occupancy and BCs
    Oxy = full(fsparse([extdof; adof],1, ...
        [ones(size(extdof)); ...
        -cons*full(U(adof)./dM(adof))], ... %change for dependence of U, original -cons*full(max(U(adof),0)./dM(adof))],
        [size(OLa.X,1) 1]));
    Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));
    
    if i==1
        Oxysave{1} = Oxy(Adof);
    end
    
    % intensities of possible events
    
    %   % (1) moving boundary DOFs
      [ii,jj_] = find(N(bdof_m,Adof)); % neighbours...
      keep = find(U(Adof(jj_)) == 0);  % ...to move to
      ikeep= ii(keep);
      ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1); %???varför??
      % remove any possibly remaining negative rates
      grad = fsparse(ii,1,max(Pr(bdof_m_(ii))-Pr(jj_),0).* ...
                     1,...%Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                     numel(bdof_m));
      moveb = full(gradquotient*grad);
    %
    % (2) also certain sources may move by the same physics
    [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
    keep = find(U(Adof(jj_)) < 1);   % ...to move to
    ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
    % remove any possibly remaining negative rates
    grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0).* ...
        1, ... %Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1)
        numel(sdof_m)); % (abs as U could be -1)
    moves = full(gradquotient*grad);
    
    %     moves=moveb;
    % (3) proliferation/death/degradation rates
%     r_prol*(U(Adof) > 0 & U(Adof) < 2).*(Oxy(Adof) > cutoff_prol);
    birth = full(r_prol*(U(Adof) > 0 & U(Adof) < 2).*(Oxy(Adof) > cutoff_prol)); %U(Adof) < 2 sätter gräns för när den får proliferate, 1?
    total_birth = sum(birth);
    birth = total_birth/total_birth * birth;
    birth(isnan(birth)) = 0;
    % (as we get some 0/0 terms if total_birth == 0);
    
    death = full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die));
    degrade = full(r_degrade*(U_dead > 0));
    %
    %   %Gillespies algorithm
    %    intens = [moveb; moves; birth; abs(death); abs(degrade)];
    intens = [birth; death; degrade; moves; moveb];
    lambda = sum(intens);
%     if i==1
%         dt=1;
%     else
%         lambda = (sum(dead_conc) + sum(prol_conc))/(length(dead_conc)+length(prol_conc));
%         lambda(isnan(lambda)) = 1;
% 
%     end
    %dt = min((1/r_degrade)*0.9, (1/r_die)*0.9);%/lambda;

    dt = 1/lambda;
%     dt=Tend/50;
    
    % report back håll some koll
    if tspan(i+1) < tt+dt
        iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
        Usave(i+1:iend) = {U};
        Udsave(i+1:iend) = {U_dead};
        Idofsave(i+1:iend) = {U(Idof)};
        
        Oxysave(i+1:iend) = {Oxy(Adof)};
        bdofsave(i+1:iend) = {bdof_m};
        sdofsave(i+1:iend) = {sdof_m};
        
        i = iend;
        
        % monitor the maximum outlier cell:
        max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));
        
        % the number of cells
        num_cells = sum(abs(U));
        
        % the rates
        %     inspect_rates = [sum(moveb) sum(moves) ...
        %                      sum(birth) sum(death) sum(degrade)];
    end
    
    %movements-----------------------------------
    D = 1; 
    %sdof_m
    rates_sdof = zeros(length(Adof),1);
    rates_test = zeros(length(Adof),1);
%     action_vec =zeros(length(Adof),1);
    
    for ind=1:length(sdof_m_)
%         Ne.moves = Ne.moves+1;
        ix = sdof_m(ind);
        ix_ = sdof_m_(ind); %index i Adof

        jx_ = find(N(ix,Adof)); %index i Adof
        Pr_diff = max(Pr(ix_)-Pr(jx_),0);
        rates_sdof(jx_) = rates_sdof(jx_) + D*Pr_diff;
        %rates_sdof(ix_) = rates_sdof(ix_) - sum(D*Pr_diff);
        rates_sdof(ix_) = rates_sdof(ix_) - sum(D*Pr_diff);
        
%         action_vec(jx_)=1;
%         action_vec(ix_)= -1;
    end
    
    %%%%%
%     jx_ = find(N(sdof_m,Adof));
%     rates_test(jx_) = D*max(Pr(sdof_m_)-Pr(jx_),0);
%     rates_test(sdof_m_)= -sum(D*max(Pr(sdof_m_)-Pr(jx_),0));
    
    
    
    %Euler for sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt; %action_vec*
    
%     U_temp = U_new(sdof_m);
%     U_temp(U_temp<cutoff) = 0;
%     U_new(sdof_m) = U_temp;
    
    %bdof
    rates_bdof = zeros(length(Adof),1);
%     action_vec =zeros(length(Adof),1);
%     Ne.moveb = Ne.moveb+1;
    for ind=1:length(bdof_m_)
        % movement of a boundary (singly occupied) voxel
        ix = bdof_m(ind);
        ix_ = bdof_m_(ind);

        jx_ = find(N(ix,Adof));
        % (will only move into an empty voxel:)
        jx_ = jx_(U(Adof(jx_)) == 0 & U_dead(Adof(jx_)) == 0);
        %rates(:,i) = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
        Pr_diff = max(Pr(ix_)-Pr(jx_),0);
        rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff;
        rates_bdof(ix_) = -sum(D*Pr_diff); 
        
%         action_vec(jx_)=1;
%         action_vec(ix_)= -1;
    end
    
    %Euler for bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt; %action_vec*
        
    %ta inte bort, höj cutoff ist
%     U_temp = U_new(bdof_m);
%     U_temp(U_temp<cutoff) = 0;
%     U_new(bdof_m) = U_temp;

    updLU = true; % boundary has changed
    
    %proliferation------------------------------
%     Ne.birth = Ne.birth+1;
%     birth_count = birth_count+1;
    

    ind_prol = find((Oxy > cutoff_prol));
    prol_conc = r_prol*U(ind_prol);
    U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;
    
    %death--------------------------------------
    % full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die))
    ind_die = find(Oxy < cutoff_die); %index for dying cells
    
%     Oxy(ind_die);
    dead_conc = r_die*U(ind_die);
    U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
    U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
    
    ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
    U_new(ind_cutoff) = 0;
    
    %Ne.death = Ne.death+1;
    
   
    %degradation--------------------------------
    %Ne.degrade = Ne.degrade+1;   
    U_deadnew(ddof) = U_dead(ddof)*(1 - r_degrade*dt); %tidssteg litet->kan inte bli negativt
    % remove degraded cells
    U_deadnew(U_dead < cutoff_deg) = 0; %ta bort för små
%     testU=U(7320)
%     testUd=U_dead(7320)
%     testOxy=Oxy(7320)
    
    updLU = true; % boundary has changed
    %end
    
    tt = tt+dt;
    report(tt,U,'');
    
    % update the visited sites
    %VU = VU | U;
end
report(tt,U,'done');

% return;

%% 

% create a GIF animation
% figure(6)
% population appearance sdof, bdof visualization
Mdof = struct('cdata',{},'colormap',{});
figure(3), clf,

Umat=cell2mat(Usave);
cmat = full(Umat/max(max(Umat)));
colorbar
caxis([0 1])
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
        
%     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 0.5 & Usave{i} <= 1);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 1 & Usave{i} <= 1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 0]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 1.5 & Usave{i} <= 2);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 1]); %[1 0 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 2 );
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 1 0]); %[1 0 Usave{i}/max(Usave{i})]
    
    patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
        'FaceColor','cyan'); %[1 0 Usave{i}/max(Usave{i})] 
    
    patch('Faces',R(sdofsave{i},:),'Vertices',V, ...
        'FaceColor','magenta'); %[1 0 Usave{i}/max(Usave{i})]      
        
%     ii = find(Usave{i} > 3 );
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 0]); %[1 0 Usave{i}/max(Usave{i})]
    
    ii = find(Usave{i} == 0 & Udsave{i} >0);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    drawnow;
    Mdof(i) = getframe(gcf);
%     pause(3)
end

% saves the GIF
movie2gif(Mdof,{Mdof([1:2 end]).cdata},'TumourMdof10.gif', ...
          'delaytime',0.5,'loopcount',0);

%%
normalfigure=1;
if normalfigure==1
% create a GIF animation

% population appearance normal
Mnormal = struct('cdata',{},'colormap',{});
figure(11), clf,

Umat=cell2mat(Usave);
cmat = full(Umat/max(max(Umat)));
colorbar
caxis([0 1])
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
    
    %     ii = find(Usave{i}>0 & Usave{i} <= 1);
    %     cvec = full(Usave{i}(ii)/max(Usave{i}(ii)));
    %     patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData', cvec);% ...
    %         %'FaceColor',[0 0 Usave{i}/max(Usave{i})]);    

    %     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 0 1]); %; [0 1 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 0.5 & Usave{i} <= 1);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 0.5 1]); %; [0 1 Usave{i}/max(Usave{i})]

%     ii = find(Usave{i} > 1);
%     cvec = full(Usave{i}(ii)/max(Usave{i}(ii)));
%     patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData', cvec);% ...
%         %'FaceColor',[0 0 Usave{i}/max(Usave{i})]);   
    
    %     ii = find(Usave{i} > 1 & Usave{i} <= 1.5);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 1 0.5]); %; [0 1 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 1.5 & Usave{i} <= 2);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0.5 1 0.2]); %[1 0 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 2 );
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[1 1 0.5]); %[1 0 Usave{i}/max(Usave{i})]

    ii = find(Usave{i} == 0 & Udsave{i} >0);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    drawnow;
    Mnormal(i) = getframe(gcf);
%     pause(2)
end
% saves the GIF
movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'TumourMnormal.gif', ...
          'delaytime',0.1,'loopcount',0);
end
%%
if deadfigure==1
% dead appearance
Mdead = struct('cdata',{},'colormap',{});
figure(10), clf,

Udmat=cell2mat(Udsave);
cmat = full(Udmat/max(max(Udmat)));
caxis([1 2])
colorbar;
%     set( h, 'YDir', 'reverse' );
colormap 'gray'

for i = 1:numel(Udsave)
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Udsave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
    
    
%     ii = find(Udsave{i} > 0 & Udsave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Udsave{i} > 0.5 & Udsave{i} <= 1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0.5 0.5 0.5]); %; [0 1 Usave{i}/max(Usave{i})] 
%   
%     ii = find(Udsave{i} >1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
    drawnow;
    Mdead(i) = getframe(gcf);
end
% saves the GIF
movie2gif(Mdead,{Mdead([1:2 end]).cdata},'animations/TumourMdead.gif', ...
          'delaytime',0.1,'loopcount',0);
end
%%

% investigate the time evolution of the different cell numbers
figure(4), clf
spsum  = @(U)(full(sum(abs(U))));
deadsum = @(U)(full(sum(U == -1)));
normsum = @(U)(full(sum(U == 1)));
prolsum = @(U)(full(sum(U == 2)));
z = cellfun(deadsum,Usave);
w = cellfun(prolsum,Usave);
x = cellfun(normsum,Usave);
y = cellfun(spsum,Usave);
p1 = plot(tspan,y);
hold on
p2 = plot(tspan,z,'k');
p3 = plot(tspan,w);
p4 = plot(tspan,x);
p3.Color = graphics_color('vermillion');
p4.Color = graphics_color('bluish green');
ylim([0 max(y)]);
xlabel('time')
ylabel('N cells')
legend('total', 'dead','double','single');

return;


%%%%%%%%%%%%%%%Extra plots

%% Plot the maxium radius through time
figure(5), clf
plot(tspan,max_radius);
xlabel('time')
ylabel('max radius')
grid on;

%% Plot the rates through time
figure(6), clf
rate_names = fieldnames(Ne);
inspect_rates_norm = inspect_rates./sum(inspect_rates,1);
bar(inspect_rates_norm','stacked','LineStyle','none') %'DisplayName',rate_names{kk});
grid on;
title('Relative and normalized rates')
xlabel('time')
ylabel('rates')
% ticks = 
set(gca, 'XTick', linspace(1,length(tspan),7))
set(gca, 'XTickLabel', round(linspace(1,tspan(end),7)))
ylim([0 1.5]);
legend(rate_names);


%% Plot Pressure
figure(7), clf,
Pr_ = full(U); Pr_(adof) = Pr(adof_);
[x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,Nvoxels));
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
surf(x_Pr_,y_Pr_,Pr_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Pr_reshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
Pr_(adof) = 0;
Pr_(idof) = Pr(idof_);
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
surf(x_Pr_,y_Pr_,Pr_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Pr_reshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)
%% Save the important data in a struct
if doSave
    saveData = struct('U', {U}, 'Usave', {Usave}, 'tspan', {tspan}, ...
        'R', {R}, 'V', {V}, 'BC1', {BC1}, 'BC2', {BC2}, ...
        'max_radius', {max_radius}, 'Ne', {Ne}, ...
        'inspect_rates', {inspect_rates}, 'alpha', {alpha}, 'Pr', {Pr}, ...
        'Adof', {Adof}, 'Nvoxels',{Nvoxels});
    filename_saveData = "saveData/saveData_T" + Tend + ...
        "_" + strjoin(string(fix(clock)),'-') + ".mat";
    save(filename_saveData, 'saveData');
end

return;

%% %% Plot U
figure(9), clf,
U_plt = full(U); U_plt(adof) = U(adof);
[x_U_plt,y_U_plt] = meshgrid(linspace(-1,1,Nvoxels));
U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
surf(x_U_plt,y_U_plt,U_pltreshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',U_pltreshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
%U_plt(adof) = 0;
U_plt(idof) = U(idof_);
U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
surf(x_U_plt,y_U_plt,U_pltreshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',U_pltreshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)

%% %% Plot Oxygen
figure(10), clf,
Oxy_ = full(U); Oxy_(adof) = Oxy(adof);
[x_Oxy_,y_Oxy_] = meshgrid(linspace(-1,1,Nvoxels));
Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Oxy_reshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
Oxy_(adof) = 0;
Oxy_(idof) = Oxy(idof_);
Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Oxy_reshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)
