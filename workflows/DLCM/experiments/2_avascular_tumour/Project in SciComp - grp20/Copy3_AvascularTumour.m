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
% for IC = [1,5]
% for alpha = [5e-2, 5e-1]
%    alpha_inv = 1/alpha; 
doGif = false;
doSave = true;
% Choose to add idof3
do_idof3 = true;
% simulation interval
% if IC == 1
    Tend = 400;
% else
%     Tend = 400;
% end
tspan = linspace(0,Tend,101);
report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

% The user specified cutoff and rate parameters for the proliferation,  
% death, degradation and consumption rules.
cons = 0.0015;        % consumption of oxygen by cells
cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
r_prol = 0.125;       % rate of proliferation of singly occupied voxels
cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
r_die = 0.125;        % rate of dea th
r_degrade = 0.01;     % rate of degradation for already dead cells

% Permeability parameters.
Drate1 = 0.01;     % into free matrix
Drate2 = 25;       % into already visited matrix 
Drate3 = 0.01;     % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];

% boundary conditions
BC1 = 10; % BC for the pressure equation for unvisited boundary
BC2 = 1; % BC for the visited boundary
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary
alpha = 1e+2;
alpha_inv = 1/alpha;

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
mesh_type = 1; % 1: cartesian, 2: hexagonal
[P,E,T,gradquotient] = basic_mesh2(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N,M] = dt_operators(P,T);
Mgamma = assemble_Mgamma(P,T);
Mgamma = Mgamma./dM;
N_Mgamma = (Mgamma ~= 0) - speye(size(Mgamma));
neigh = full(sum(N,2));

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

% Initial population
IC = 5; % Choose initial condition (1,2,3,4,5,6)
R1 = 0.35; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)
U = setInitialCondition(IC,R1,R2,P,Nvoxels);

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
Ne = struct('moveb1',0,'moveb2',0,'moves',0,'birth',0,'death',0,'degrade',0);
timing_vec = zeros(6,length(tspan)+1);

% max radius and inspect_rates vectors
max_radius = zeros(1,numel(tspan));
inspect_rates = zeros(length(fieldnames(Ne)),numel(tspan));

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

% tic
while tt <= tspan(end)
  %% classify the DOFs
  adof = find(U); % all filled voxels
  % singularly occupied voxels on the boundary:
  bdof_m = find(N*(U ~= 0) < neigh & abs(U) == 1);
  sdof = find(U > 1); % voxels with 2 cells
  % voxels with 2 cells in them _which may move_, with a voxel
  % containing less number of cells next to it (actually 1 or 0):
  sdof_m = find(N*(U > 1) < neigh & U > 1);
  Idof = (N*(U ~= 0) > 0 & U == 0); % empty voxels touching occupied ones
  idof1 = find(Idof & ~VU); % "external" OBC1
  idof2 = find(Idof & VU);  % "internal" OBC2
  idof3 = find(~VU & N_Mgamma*VU > 0); % boundary around "visited voxels"
  idof3 = setdiff(idof3,idof1);
  idof = find(Idof);
  % Divide bdof_m into the singularly occupied voxels
  % who touch either boundary voxels from idof1 or idof2 
  [ind_r,ind_c] = find(N(bdof_m,:));
  bdof_m1 = bdof_m(unique(ind_r(find(ismember(ind_c,idof1)))));
  bdof_m2 = bdof_m(unique(ind_r(find(ismember(ind_c,idof2)))));
  bdof_m1(find(ismember(bdof_m1, bdof_m2))) = [];
  bdof_m1 = reshape(bdof_m1,[],1);
  bdof_m2 = reshape(bdof_m2,[],1);

  if do_idof3
    idof1 = [idof1;idof3];
    idof = [idof;idof3];
  end

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];

  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';  
 [bdof_m_,bdof_m1_,bdof_m2_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_] = ...
      map(Adof_,Adof,bdof_m,bdof_m1,bdof_m2,sdof,sdof_m,idof1,idof2,idof,adof);
  %% Update LU
  if updLU
     % pressure Laplacian
    La.X = L(Adof,Adof);

    %%% ADD BC to LHS
    % Robin to external idof (idof1)
    Lai1 = fsparse(idof1_,idof1_,1,size(La.X));
    a_Lai1 = speye(size(Lai1)) - Lai1;
    
    La.X = La.X - Lai1*La.X + Lai1;

    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
    updLU = false; % assume we can reuse
  end

  %% Caculate laplacians

  % RHS source term proportional to the over-occupancy and BCs
  Pr = full(fsparse(sdof_,1,1./dM(sdof), ...
                  [size(La.X,1) 1]));     % RHS first...

  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution
%   Pr = normalize(Pr,'range');

  % RHS source term proportional to the over-occupancy and BCs
  Oxy = full(fsparse([extdof; adof],1, ...
                     [ones(size(extdof)); ... 
                      -cons*full(max(U(adof),0)./dM(adof))], ...
                     [size(OLa.X,1) 1]));
  Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));

  %% Measure events probabilites
  % intensities of possible events

  % (1.1) moving boundary DOFs to idof1
  [ii,jj_] = find(N(bdof_m1,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) == 0);  % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(bdof_m1_(ii))-Pr(jj_),0).* ...
                 Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                 numel(bdof_m1));
  moveb1 = full(gradquotient*grad);

  % (1.2) moving boundary DOFs to idof2
  [ii,jj_] = find(N(bdof_m2,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) == 0);  % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(bdof_m2_(ii))-Pr(jj_),0).* ...
                 Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                 numel(bdof_m2));
  moveb2 = full(gradquotient*grad);

  % (2) also certain sources may move by the same physics
  [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) < 2);   % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0).* ...
               Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1), ...
                 numel(sdof_m)); % (abs as U could be -1)
  moves = full(gradquotient*grad);

  % (3) proliferation/death/degradation rates
  birth = full(r_prol*(U(Adof) == 1).*(Oxy(Adof) > cutoff_prol));
  total_birth = sum(birth);
%   birth = total_birth/total_birth * birth;
  birth(isnan(birth)) = 0;
  % (as we get some 0/0 terms if total_birth == 0);

  death = full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die));
  degrade = full(r_degrade*(U(Adof) == -1));

  %% Caclutate which is suppose to happen
  intens = [moveb1; moveb2; moves; birth; death; degrade];
  lambda = sum(intens);
%   fprintf('sum(moveb2)/sum(moveb1) = %d \n', sum(moveb2)/sum(moveb1));
%   fprintf('(sum(moveb2) + sum(moveb1))/sum(moves) = %d \n', (sum(moveb2) + sum(moveb1))/sum(moves));
  dt = -reallog(rand)/lambda; 
  rnd = rand*lambda;
  cum = intens(1);
  ix_ = 1;
  while rnd > cum
    ix_ = ix_+1;
    cum = cum+intens(ix_);
  end
  % (now ix_ points to the intensity which fired first)

  %% Execute the event that happens
  if ix_ <= numel(moveb1)
    Ne.moveb1 = Ne.moveb1+1;
    % movement of a boundary (singly occupied) voxel
    ix_ = bdof_m1_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    % (will only move into an empty voxel:)
    jx_ = jx_(U(Adof(jx_)) == 0);
    rates = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    U(n) = U(ix);
    U(ix) = 0;
    updLU = true; % boundary has changed
  elseif ix_ <= numel(moveb1)+numel(moveb2)
    Ne.moveb2 = Ne.moveb2+1;
    % movement of a boundary (singly occupied) voxel
    ix_ = ix_-numel(moveb1);
    ix_ = bdof_m2_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    % (will only move into an empty voxel:)
    jx_ = jx_(U(Adof(jx_)) == 0);
    rates = Drate_(2*VU(Adof(jx_))+1).*max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    U(n) = U(ix);
    U(ix) = 0;
    updLU = true; % boundary has changed
  elseif ix_ <= numel(moveb1)+numel(moveb2)+numel(moves)
    Ne.moves = Ne.moves+1;
    % movement of a cell in a doubly occupied voxel
    ix_ = ix_-numel(moveb1)-numel(moveb2);
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
  elseif ix_ <= numel(moveb1)+numel(moveb2)+numel(moves)+numel(birth)
    Ne.birth = Ne.birth+1;
    % proliferation
    birth_count = birth_count+1;
    ix_ = ix_-numel(moveb1)-numel(moveb2)-numel(moves);
    ix = Adof(ix_);
    U(ix) = U(ix)+1;
  elseif ix_ <= numel(moveb1)+numel(moveb2)+numel(moves)+numel(birth)+numel(death)
    Ne.death = Ne.death+1;
    % death
    ix_ = ix_-numel(moveb1)-numel(moveb2)-numel(moves)-numel(birth);
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
    ix_ = ix_-numel(moveb1)-numel(moveb2)-numel(moves)-numel(birth)-numel(death);
    ix = Adof(ix_);
    U(ix) = 0;
    updLU = true; % boundary has changed
  end

  %% Rest of while loop
  % report back
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};

    % monitor the maximum outlier cell:
    max_radius(i+1:iend) = sqrt(max(P(1,adof).^2+P(2,adof).^2));


    % the number of cells
    num_cells = sum(abs(U));

    % the rates
    inspect_rates(:,i) = [sum(moveb1) sum(moveb2) sum(moves) ...
                     sum(birth) sum(death) sum(degrade)];

    i = iend;
    
%     figure(8);
%     plot(Pr)
%     hold on;
%     xxx = 1:length(Pr);
%     plot(1:length(Pr),ones(size(xxx))*mean(Pr), 'g', 'LineWidth',1.5)
%     legend('Pr', sprintf('Mean Pr = %0.3f',mean(Pr)));
  end

  tt = tt+dt;
  report(tt,U,'');

  % update the visited sites
  VU = VU | U;
end
% toc
report(tt,U,'done');

% return;
%%
% create a GIF animation

% population appearance
figure(3), clf,
if doGif
    M = struct('cdata',{},'colormap',{});
    for i = 1:2:numel(Usave)
      patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
      hold on,
      axis([-1 1 -1 1]); axis square, axis off
      ii = find(Usave{i} == 1);
      patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('bluish green'));
      ii = find(Usave{i} == 2);
      patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('vermillion'));
      ii = find(Usave{i} == -1);
      patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
      title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
      drawnow;
      M(i) = getframe(gcf);
    end
else
    patchCurrentCells;
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


%% Plot the maxium radius through time
% figure(5), clf
% plot(tspan,max_radius);
% xlabel('time')
% ylabel('max radius')
% grid on;

%% Plot the rates through time
figure(6), clf
plotRates;

%% Plot Pressure
figure(7), clf,
if 1
    plotPressure;
else
    % Here is some weird try on a 3D bar plot
    % OBSERVE: takes a long time to plot
    plotPressureBars;
end
    
%% Save the important data in a struct
if doSave
    saveData = struct('U', {U}, 'VU', {VU}, 'Usave', {Usave}, 'tspan', {tspan}, ...
        'R', {R}, 'V', {V}, 'N', {N}, 'BC1', {BC1}, 'BC2', {BC2}, ...
        'max_radius', {max_radius}, 'Ne', {Ne}, 'inspect_rates', {inspect_rates}, ...
        'alpha', {alpha}, 'Pr', {Pr}, 'Adof', {Adof},  ...
        'adof', {adof}, 'adof_', {adof_}, 'idof', {idof}, 'idof_', {idof_}, ...
        'Nvoxels',{Nvoxels}, 'IC', {IC}, 'R1', {R1}, 'R2', {R2},...
        'P', {P}, 'bdof_m', {bdof_m}, 'bdof_m_', {bdof_m_}, ...
        'sdof_m', {sdof_m},'sdof_m_', {sdof_m_}, 'gradquotient', {gradquotient}, ...
        'Tend', {Tend}, 'mesh_type', {mesh_type}, 'Drate_', {Drate_}, 'do_idof3', {do_idof3});
    filename_ = "alpha" + erase(sprintf('%0.0e',alpha),'.');
    if mesh_type == 2
        filename_ = filename_ + "_HEX";
    end
    filename_ = filename_ + "_" + strjoin(string(fix(clock)),'-');
    filename_saveData = "saveData/saveData_" + filename_ + ".mat";
    save(filename_saveData,'-struct','saveData');
    % print('-f7', "images/pressurePlot_" + filename_, '-depsc');
end
% end
% end
return;

% % saves the GIF
% movie2gif(M,{M([1:2 end]).cdata},'animations/Tumour.gif', ...
%           'delaytime',0.1,'loopcount',0);
