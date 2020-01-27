% Basic test of simulation of a population of cells, 3D version.

% S. Engblom 2017-08-29 (revision, cleanup)
% S. Engblom 2016-11-13 (revision, 3D version)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

% cells live in a square of Nvoxels-by-Nvoxels-by-Nvoxels
Nvoxels = 20;

% Cartesian (midpoint-)discretization
h = 2/Nvoxels;
X = -1+h/2:h:1-h/2;
[X,Y,Z] = ndgrid(X,X,X);
X = X(:);
Y = Y(:);
Z = Z(:);
P = [X Y Z]';

% assemble minus the Laplacian on this grid (ignoring BCs)
L1 = [0 0 0; 0 1 0; 0 0 0];
L2 = [0 1 0; 1 -6 1; 0 1 0];
L = cat(3,L1,L2,L1)/h^2;
L = -ndop(L,[2 2 2],[Nvoxels Nvoxels Nvoxels]);
dM = h^3;

% neighbor matrix
N1 = [0 0 0; 0 1 0; 0 0 0];
N2 = [0 1 0; 1 0 1; 0 1 0];
N = cat(3,N1,N2,N1);
N = ndop(N,[2 2 2],[Nvoxels Nvoxels Nvoxels]);

% initial populations of cells
ii = find(max(abs(P),[],1) <= 0.4);
U = fsparse(ii(:),1,2,[Nvoxels^3 1]);
% (format: U = {0,1,2} for {empty,one cell,two cells})

% simulation interval
Tend = 20;
% solution recorded at this grid:
tspan = linspace(0,Tend,41);

report(tspan,U,'init');

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));

tt = tspan(1);
Usave{1} = U;
i = 1;
while tt <= tspan(end)
  % classify the Degrees-Of-Freedoms (DOFs for short)

  % active DOFs: occupied voxels
  adof = find(U);

  % boundary DOFs _which may move_: containing one cell per voxel and
  % with an empty voxel besides it
  bdof_m = find(N*(U > 0) < 6 & U == 1);

  % source DOFs containing more than one cell per voxel (actually 2)
  sdof = find(U > 1);

  % source DOFs _which may move_: with a voxel containing less number of
  % cells next to it (actually 1 or 0)
  sdof_m = find(N*(U > 1) < 6 & U > 1);

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

  % the currently active Laplacian
  La = L(Adof,Adof);
  % selecting all injection DOFs
  Lai = fsparse(idof_,idof_,1,size(La));
  % remove eqs (rows) for injection DOFs, replace with direct injection
  La = La-Lai*La+Lai;

  % RHS source term proportional to the over-occupancy (the rest of the
  % zeros ensure the homogeneous Dirichlet BC is satisfied at idof)
  Pr = La\full(fsparse(sdof_,1,1/dM,[size(La,1) 1]));

  % Movement rates are proportional to minus the gradient (also with
  % homogeneous Dirichlet boundary conditions in the layer of voxels
  % just outside the boundary).

  % compute intensities of possible events

  % moving boundary DOFs: each empty neighbour voxel is associated with
  % a flow rate proportional to the pressure gradient in that
  % direction integrated over the corresponding edge
  grad = sum(N(bdof_m,idof),2).*Pr(bdof_m_);  % (note: Dirichlet 0 BCs)
  moveb = full(grad);

  % also certain sources may move by the same physics
  [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) < 2);   % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0),numel(sdof_m));
  moves = full(grad);

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
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
  end
  tt = tt+dt;
  report(tt,U,'');
end
report(tt,U,'done');

return;

% a few frames
j = 1;
for i = [1 2 numel(Usave)]
  j = j+1;
  figure(j), clf, hold on,
  ii = find(Usave{i} == 1);
  keep = find(P(1,ii)+P(2,ii)+P(3,ii) <= 0);
  ii = ii(keep);
  h = plot3(P(1,ii),P(2,ii),P(3,ii),'go');
  set(h,'MarkerEdgeColor',graphics_color('bluish green'));
  set(h,'MarkerFaceColor',graphics_color('bluish green'));
  set(h,'MarkerSize',2);
  jj = find(Usave{i} > 1);
  keep = find(P(1,jj)+P(2,jj)+P(3,jj) <= 0);
  jj = jj(keep);
  h = plot3(P(1,jj),P(2,jj),P(3,jj),'ro');
  set(h,'MarkerEdgeColor',graphics_color('vermillion'));
  set(h,'MarkerFaceColor',graphics_color('vermillion'));
  set(h,'MarkerSize',3);
  axis([-1 1 -1 1 -1 1]); axis square, 
  set(gca,'xticklabel',[])
  set(gca,'yticklabel',[])
  set(gca,'zticklabel',[])
  view(30,30);

  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 340 220]);
  drawnow;
end

% uncomment to save:
% $$$ figure(2),
% $$$ print -depsc figures/basic_test_1_3D.eps
% $$$ figure(3),
% $$$ print -depsc figures/basic_test_2_3D.eps
% $$$ figure(4),
% $$$ print -depsc figures/basic_test_end_3D.eps