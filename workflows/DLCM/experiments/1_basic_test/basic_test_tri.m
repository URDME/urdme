% Basic test of simulation of a population of cells.
%
%   Case of unstructrued triangular mesh.

% S. Engblo 2017-12-26 (revision)
% S. Engblom 2017-01-14 (_tri)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-25 (seriously finalized the physics)
% S. Engblom 2016-12-20 (finalized the physics)
% S. Engblom 2016-12-09 (hexagonal mesh)
% S. Engblom 2016-12-08 (notes on thinning)
% S. Engblom 2016-12-02 (revision)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

% diffusive pressure rate
Drate = 1;

% if geom == 3
if ~exist('reuse','var') || ~reuse
  % cells live in a circular geometry
  % petri_geom; % --> export geometry and boundary into [G,B]
  load experiments/petri_mesh/petri_geom % load this [G,B]

  % create mesh
  [P,E,T] = initmesh(G,'hmax',0.0577); % Nvoxels = 1602 ~ 40*40

  % for visualization using Voronoi patches
  [V,R] = mesh2dual(P,E,T,'voronoi');

  % each node in the mesh is a (midpoint of a) voxel
  Nvoxels = size(P,2);

  % assemble operators
  [L,dM,N] = dt_operators(P,T);

  % integration of pressure gradients over edges

  % distance between two voxel midpoints: D1(k) = norm(P(:,i(k)),P(:,j(k)))
  [i,j] = find(N);  
  D1 = sqrt(sum((P(:,i)-P(:,j)).^2))';

  % edge length: L1 = edge length of edge between voxel i and j
  L1 = zeros(size(D1));
  for k = 1:numel(i)
    n = fsetop('intersect',R(i(k),:),R(j(k),:));
    n(isnan(n)) = [];
    L1(k) = norm(diff(V(n,:)));
  end

  % quotient (this is not a scalar here!)
  gradquotient = fsparse(i,j,L1./D1,size(L));
end

% initial population
ii = find(max(abs(P),[],1) <= 0.4);
U = fsparse(ii(:),1,2,[Nvoxels 1]);
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
neigh = full(sum(N,2));
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
  [ii,jj_,gq] = find(gradquotient(bdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) == 0);                % ...to move to
  ii = reshape(ii(keep),[],1);
  gq = reshape(gq(keep),[],1);
  grad = fsparse(ii,1,gq.*Pr(bdof_m_(ii)),numel(bdof_m));
  moveb = Drate*full(grad);

  % also certain sources may move by the same physics
  [ii,jj_,gq] = find(gradquotient(sdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) < 2);                 % ...to move to
  ii = reshape(ii(keep),[],1);
  jj_ = reshape(jj_(keep),[],1);
  gq = reshape(gq(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,gq.*max(Pr(sdof_m_(ii))-Pr(jj_),0),numel(sdof_m));
  moves = Drate*full(grad);

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
    ix_ = bdof_m_(ix_);
    ix = Adof(ix_);

    % select neighbor to move to
    jx_ = find(N(ix,Adof));
    jx_ = jx_(U(Adof(jx_)) == 0); % only allow moves to less populated voxels
    % the pressure gradient drives the movement
    rates = full(gradquotient(ix,Adof(jx_)))'.*Pr(ix_);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));
    
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
    rates = full(gradquotient(ix,Adof(jx_)))'.*max(Pr(ix_)-Pr(jx_),0);
    % (thinning of negative gradients)
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

% last frame
if report_progress
  figure(1), clf,
  for i = numel(Usave)
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    ii = find(Usave{i} == 1);
    patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
    jj = find(Usave{i} > 1);
    patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
    drawnow;
  end
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
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
  drawnow;
  title(sprintf('t = %f',tspan(i)));
  M(i) = getframe(gcf);
end

% $$$ movie2gif(M,{M(1:2).cdata},'animations/basic_test_tri.gif', ...
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
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);

  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 340 220]);
  drawnow;  
end

% uncomment to save:
% $$$ figure(2),
% $$$ print -depsc figures/basic_test_1_tri.eps
% $$$ figure(3),
% $$$ print -depsc figures/basic_test_2_tri.eps
% $$$ figure(4),
% $$$ print -depsc figures/basic_test_end_tri.eps
