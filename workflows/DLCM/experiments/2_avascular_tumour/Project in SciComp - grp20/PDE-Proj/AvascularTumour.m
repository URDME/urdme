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

% simulation interval
Tend =50;
tspan = linspace(0,Tend,101);
% report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

%     The user specified cutoff and rate parameters for the proliferation,
%     death, degradation and consumption rules.

rates_type = 2;
if rates_type == 1
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
else
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.45;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.35;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
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
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < 0.09); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
else
    % initial population: circular blob of dead cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < 0.05); % radius of the initial blob
    U = fsparse(ii(:),1,0,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,1,[Nvoxels^2 1]);
end

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;
Udsave = cell(1,numel(tspan));
Udsave{1} = U_dead;

Oxysave = cell(1,numel(tspan));


birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   %add Dirichlet(?) oxygen at the oxygen source/outer circle
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

while tt <= tspan(end)
    %   Ncells=full(sum(Usave{i}~=0))
    %   Ndead=full(sum(Usave{i}==-1))
    %   Nalive=full(sum(Usave{i}>0))
    
    % classify the DOFs
    adof = find(U|U_dead); % all filled voxels (U_dead not necessary if U=-1)
    % singularly occupied voxels on the boundary: (will not be a source of
    % moving cells)?
    bdof_m = find(N*(U ~= 0) < neigh & (U > 0 & U <= 1));
    sdof = find(U > 1); % voxels with 2 cells
    % voxels with 2 cells in them _which may move_, with a voxel
    % containing less number of cells next to it (actually 1 or 0):
    sdof_m = find(N*(U > 1 | U_dead >1)<neigh & U > 1); %kanske större än 1
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
    [bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_] = ...
        map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof,adof);
    
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
    %   [ii,jj_] = find(N(bdof_m,Adof)); % neighbours...
    %   keep = find(U(Adof(jj_)) == 0);  % ...to move to
    %   ikeep= ii(keep);
    %   ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1); %???varför??
    %   % remove any possibly remaining negative rates
    %   grad = fsparse(ii,1,max(Pr(bdof_m_(ii))-Pr(jj_),0).* ...
    %                  Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
    %                  numel(bdof_m));
    %   moveb = full(gradquotient*grad);
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
    intens = [birth; death; degrade; moves];
    lambda = sum(intens);
    %     if i==1
    %         dt=1;
    %     else
    %         lambda = (sum(dead_conc) + sum(prol_conc))/(length(dead_conc)+length(prol_conc));
    %         lambda(isnan(lambda)) = 1;
    %
    %     end
    %     dt = min((1/r_degrade)*0.9, (1/r_die)*0.9);%/lambda;
    
    dt = 1/lambda;
    %dt=1;
    
    % report back håll some koll
    if tspan(i+1) < tt+dt
        iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
        Usave(i+1:iend) = {U};
        Udsave(i+1:iend) = {U_dead};
        
        Oxysave(i+1:iend) = {Oxy(Adof)};
        
        i = iend;
        
        % monitor the maximum outlier cell:
        max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));
        
        % the number of cells
        num_cells = sum(abs(U));
        
        % the rates
        %     inspect_rates = [sum(moveb) sum(moves) ...
        %                      sum(birth) sum(death) sum(degrade)];
    end
    
    rates = zeros(Adof,bdof_m_); 
    for i=1:length(bdof_m_)
        Ne.moveb = Ne.moveb+1;
        % movement of a boundary (singly occupied) voxel
        ix = Adof(i);

        jx_ = find(N(ix,Adof));
        % (will only move into an empty voxel:)
        jx_ = jx_(U(Adof(jx_)) == 0);
        rates(:,i) = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
%         n = Adof(jx_(m));
    end
    % execute event: move from ix to n
    U(n) = U(ix);
    U(ix) = 0;
    updLU = true; % boundary has changed

    %elseif ix_ <= numel(moveb)+numel(moves)
    %     Ne.moves = Ne.moves+1;
    %     % movement of a cell in a doubly occupied voxel
    %     ix_ = ix_-numel(moveb);
    %     ix_ = sdof_m_(ix_);
    %     ix = Adof(ix_);
    %
    %     jx_ = find(N(ix,Adof));
    %     % (won't move into a voxel containing a dead -1 cell:)
    %     jx_ = jx_(-1 < U(Adof(jx_)) & U(Adof(jx_)) < 2);
    %     rates = Drate_(2*VU(Adof(jx_))+abs(U(Adof(jx_)))+1).* ...
    %             max(Pr(ix_)-Pr(jx_),0);
    %     m = find(cumsum(rates) > rand*sum(rates),1,'first');
    %     n = Adof(jx_(m));
    
    %       N_sdof = N;
    %       N_sdof(sdof_m,:)=2;
    %       s_mat = N_sdof-N;
    %
    %       %s_mat(sdof_m,:) = 0;
    %       Pr_mat = s_mat(Adof,Adof).*Pr - s_mat(Adof,Adof).*Pr';
    %       Pr_diff = sum(Pr_mat,2);
    D = 1; %Fix later, permeability etc D_rate
    U(Adof) = U(Adof) - D*grad*dt;
    %     % execute event: move from ix to n
    %     if U(n) == 0, updLU = true; end % boundary has changed
    %     U(n) = U(n)+1;
    %     U(ix) = U(ix)-1;
    
    
    %   elseif ix_ <= numel(moveb)+numel(moves)+numel(birth)
    
    %Ne.birth = Ne.birth+1;
    % proliferation
    %birth_count = birth_count+1;
    
    ind_prol = find((Oxy > cutoff_prol));
    prol_conc = r_prol*U(ind_prol);
    U(ind_prol)=U(ind_prol)+prol_conc*dt;
    %     ind_test =  find(U>4.1);
    %     U(ind_test) = 4;
    
    
    %death--------------------------------------
    % full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die))
    ind_die = find(Oxy < cutoff_die); %index for dying cells
    
    Oxy(ind_die);
    dead_conc = r_die*U(ind_die);
    U(ind_die) = U(ind_die) - dead_conc*dt;
    U_dead(ind_die) = U_dead(ind_die) + dead_conc*dt;
    
    ind_cutoff =  find(U<0.2 & (Oxy < cutoff_die));
    U(ind_cutoff) = 0;
    
    %Ne.death = Ne.death+1;
    
    %degradation--------------------------------
    %Ne.degrade = Ne.degrade+1;
    
    
    U_dead(ddof) = U_dead(ddof)*(1 - r_degrade*dt); %tidssteg litet->kan inte bli negativt
    % remove degraded cells
    U_dead(U_dead < 0.0001) = 0; %ta bort för små
    testU=U(7320)
    testUd=U_dead(7320)
    testOxy=Oxy(7320)
    
    
    updLU = true; % boundary has changed
    %end
    
    tt = tt+dt;
    report(tt,U,'');
    
    % update the visited sites
    %VU = VU | U;
end
report(tt,U,'done');

% return;

% create a GIF animation
% figure(6)
% population appearance
M = struct('cdata',{},'colormap',{});
figure(3), clf,
% for i = 1:numel(Usave)
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    ii = find(Usave{i} > 0 & Usave{i} <= 1);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green')); %; [0 1 Usave{i}/max(Usave{i})]
    ii = find(Usave{i} > 1);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion')); %[1 0 Usave{i}/max(Usave{i})]
    ii = find(Usave{i} == 0 & Udsave{i} >0);
    %   test=Udsave{i}(7320)
    %   color = full(Udsave{i}(ii))/max(Udsave{i}(ii));
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
    drawnow;
    M(i) = getframe(gcf);
end

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

% saves the GIF
% movie2gif(M,{M([1:2 end]).cdata},'animations/Tumour.gif', ...
%           'delaytime',0.1,'loopcount',0);
