% Basic test of simulation of a population of cells.
%
%   Relaxation to equilibrium: here the cells start distributed within
%   an internal square region filled with 2 cells (red) per voxel. The
%   system then relaxes to equilibrium by, at any given point in time
%   assume quasi steady-state and solve an equation for the 'cellular
%   pressure'.

%   Outline of method: at any given point in time we assume quasi
%   steady-state and solve an equation for the 'cellular pressure' in
%   the form of a Laplacian with source terms. Cells can move only
%     -when they have empty voxels next to them (i.e. at boundary
%     points), and here the gradient of the pressure in that direction
%     is understood as a rate per unit of time to change position,
%     -or in general, when a neighbor voxel is less populated than the
%     current one, then the (positive) gradient of the pressure in
%     that direction is again understood as a rate per unit of time to
%     move.
%
%   This stochastic process is simulated in continuous time in the
%   form of a Markov chain.
%
%   The algorithm thus consists of two steps:
%     (1) assume quasi equlibrium (no cells make any large movements,
%     only small movements about each center of mass) - solve the
%     pressure equation with sources at each cell position where there
%     are > 1 cells,
%     (2) all rates determined in this way now imply a cell which can
%     move - find out which ones moves first, and move it.

% S. Engblom 2017-12-20 (revision, more cleanup)
% S. Engblom 2017-08-29 (revision, cleanup)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-25 (seriously finalized the physics)
% S. Engblom 2016-12-20 (finalized the physics)
% S. Engblom 2016-12-09 (hexagonal mesh)
% S. Engblom 2016-12-08 (notes on thinning)
% S. Engblom 2016-12-02 (revision)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 40;
% select Nvoxels large enough so that no cell touches the boundary

% diffusive pressure rate
Drate = 1;

% fetch discretization (mesh_type = 1 or 2)
if ~exist('mesh_type','var'), error('Must define mesh_type.'); end
[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);

% initial population
ii = find(max(abs(P),[],1) <= 0.4);
U = fsparse(ii(:),1,2,[Nvoxels^2 1]);
% (format: U = {0,1,2} for {empty,one cell,two cells})

% simulation interval
Tend = 100;
% solution recorded at this grid:
tspan = linspace(0,Tend,41);

if ~exist('report_progress','var')
  report_progress = true;
end
if report_progress
  report(tspan,U,'init');
else
  report(tspan,U,'none');
end

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));

tt = tspan(1);
Usave{1} = U;
i = 1;
neigh = full(max(sum(N,2))); % (= 4 or 6)
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
while tt <= tspan(end)
% $$$   % visualization (somewhat slow)
% $$$   figure(1), clf,
% $$$   patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
% $$$         'EdgeColor','none');
% $$$   hold on,
% $$$   axis([-1 1 -1 1]); axis square, axis off
% $$$   ii = find(U == 1);
% $$$   patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
% $$$   jj = find(U > 1);
% $$$   patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
% $$$   drawnow;

  % classify the Degrees-Of-Freedoms (DOFs for short)

  % active DOFs: occupied voxels
  adof = find(U);

  % boundary DOFs _which may move_: containing one cell per voxel and
  % with an empty voxel besides it
  bdof_m = find(N*(U > 0) < neigh & U == 1);

  % source DOFs containing more than one cell per voxel (actually 2)
  sdof = find(U > 1);

  % source DOFs _which may move_: with a voxel containing less number of
  % cells next to it (actually 1 or 0)
  sdof_m = find(N*(U > 1) < neigh & U > 1);

  % injection DOFs: the layer of voxels "just outside" adof where we may
  % inject a BC (pressure 0 of the free matrix)
  idof = find(N*(U ~= 0) > 0 & U == 0);

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];
  % (after this adof is not used)
  
  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';  
  [bdof_m_,sdof_,sdof_m_,idof_] = ...
      map(Adof_,Adof,bdof_m,sdof,sdof_m,idof);
  
  if updLU
    % pressure Laplacian - factorization is reused
    La.X = L(Adof,Adof);

    % selecting all injection DOFs
    Lai = fsparse(idof_,idof_,1,size(La.X));
    % remove eqs (rows) for injection DOFs, replace with direct injection
    La.X = La.X-Lai*La.X+Lai;

    % factorize
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
    updLU = false; % assume we can reuse
  end

  % RHS source term proportional to the over-occupancy (the rest of the
  % zeros ensure the homogeneous Dirichlet BC is satisfied at idof)
  Pr = full(fsparse(sdof_,1,1./dM(sdof),[size(La.X,1) 1])); % RHS
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr));
  % with [L,U,p,q,R] = lu(A,'vector'), then R(:,p)\A(:,q) = L*U  
  %Pr = La.X\full(fsparse(sdof_,1,1./dM(sdof),[size(La.X,1) 1]));

  % Movement rates are proportional to minus the pressure gradient (with
  % homogeneous Dirichlet boundary conditions in the layer of voxels
  % just outside the boundary).

  % compute intensities of possible events

  % moving boundary DOFs: each empty neighbour voxel is associated with
  % a flow rate proportional to the pressure gradient in that
  % direction integrated over the corresponding edge
  grad = sum(N(bdof_m,idof),2).*Pr(bdof_m_);  % (note: Dirichlet 0 BCs)
  moveb = Drate*full(gradquotient*grad);

  % also certain sources may move by the same physics
  [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) < 2);   % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0),numel(sdof_m));
  moves = Drate*full(gradquotient*grad);

  intens = [moveb; moves];

  % waiting time until next event and the event itself
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

  % record values as needed
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    i = iend;
  end

  if ix_ <= numel(moveb)
    % movement of a boundary (singly occupied) voxel: find DOF of event in
    % global enumeration
    ix = Adof(bdof_m_(ix_));

    % neighbors of boundary DOF ix... which are empty... which is the one
    % selected at random
    n = find(N(ix,:));
    n = n(U(n) == 0);
    n = n(ceil(numel(n)*rand));

    % execute event: move from ix to n
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
    updLU = true;
  else
    % movement of a cell in a doubly occupied voxel
    ix_ = ix_-numel(moveb);
    ix_ = sdof_m_(ix_);
    ix = Adof(ix_);

    % select neighbor to move to
    jx_ = find(N(ix,Adof));
    jx_ = jx_(U(Adof(jx_)) < 2); % only allow moves to less populated voxels
    % the pressure gradient drives the movement
    rates = max(Pr(ix_)-Pr(jx_),0); % thinning of negative gradients
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n
    if U(n) == 0
      updLU = true;
    end
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
  end
  tt = tt+dt;
  report(tt,U,'');
end
report(tt,U,'done');

% P2A measure of circularity over time
P2A = zeros(1,numel(Usave));
for i = 1:numel(Usave)
  ii = find(Usave{i});
  [j,Area] = convhull(P(1,ii),P(2,ii));
  Perimeter = sum(sqrt(diff(P(1,ii(j))).^2+diff(P(2,ii(j))).^2));
  P2A(i) = Perimeter^2/(4*pi*Area);
end

return;

% postprocessing: visualization in true time
M = struct('cdata',{},'colormap',{});
figure(1),
for i = 1:numel(Usave)
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  axis([-1 1 -1 1]); axis square, axis off
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'));
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'));
  drawnow;
  title(sprintf('t = %f',tspan(i)));
  M(i) = getframe(gcf);
end

% P2A measure of circularity over time
P2A = zeros(1,numel(Usave));
for i = 1:numel(Usave)
  ii = find(Usave{i});
  [j,Area] = convhull(P(1,ii),P(2,ii));
  Perimeter = sum(sqrt(diff(P(1,ii(j))).^2+diff(P(2,ii(j))).^2));
  P2A(i) = Perimeter^2/(4*pi*Area);
end
figure(2), plot(tspan,P2A);
xlabel('t');
ylabel('P2A');

% uncomment to save:
% $$$ movie2gif(M,{M(1:2).cdata},'animations/basic_test.gif', ...
% $$$           'delaytime',0.1,'loopcount',0);
% $$$ movie2gif(M,{M(1:2).cdata},'animations/basic_test_hex.gif', ...
% $$$           'delaytime',0.1,'loopcount',0);

% a few frames
j = 1;
for i = [1 2 numel(Usave)]
  j = j+1;
  figure(j), clf,
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  axis([-1 1 -1 1]); axis square, axis off
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'));
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'));

  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 340 220]);
  drawnow;  
end

% uncomment to save:
% $$$ figure(2),
% $$$ print -depsc figures/basic_test_1.eps
% $$$ figure(3),
% $$$ print -depsc figures/basic_test_2.eps
% $$$ figure(4),
% $$$ print -depsc figures/basic_test_end.eps

% $$$ figure(2),
% $$$ print -depsc figures/basic_test_1_hex.eps
% $$$ figure(3),
% $$$ print -depsc figures/basic_test_2_hex.eps
% $$$ figure(4),
% $$$ print -depsc figures/basic_test_end_hex.eps
