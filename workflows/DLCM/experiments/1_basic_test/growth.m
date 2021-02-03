% Basic growth model.

% S. Engblom 2017-12-26 (revision)
% S. Engblom 2017-03-30 (adapted from delta_notch)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-27 (finalized physics)
% S. Engblom 2016-12-21 (updated the physics, changed to hex-cells)
% S. Engblom 2016-11-10

rng(222); % repeatable results

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 40;

% diffusive pressure rate
Drate = 2;

% build the hex-mesh
[P,E,T,gradquotient] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);

% birth process
lam = 1;    % (initial) proliferation rate
rho = 0.10; % radius of proliferation

% initial population
ii = find(sqrt(sum(P.^2)) <= 2*rho);
U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
% sparse solution vector U(i) = {0,1,2} cells for voxel i

% *** not implemented "all the way" - currently unused:
% identification index
ID = fsparse(ii(:),1,(1:numel(ii))',[Nvoxels^2 2]);
% Format: this is a sparse solution vector over the full mesh [P,E,T],
% where the second column is unused until cells proliferate, a filled
% row then means the corresponding voxel contains two cells (which is
% the maximum). ID is currently unused but could be employed to keep
% track of individual cells.

% simulation interval
Tend = 10;
% solution recorded at this grid:
tspan = linspace(0,Tend,41);

report(tspan,U,'init');

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
IDsave = cell(1,numel(tspan));

tt = tspan(1);
Usave{1} = U;
IDsave{1} = ID;
i = 1;
neigh = 6;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
while tt <= tspan(end)
% $$$   % visualization (somewhat slow)
% $$$   figure(1), clf,
% $$$   patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9]);
% $$$   hold on,
% $$$   axis([-1 1 -1 1]); axis square, axis off
% $$$   ii = find(U == 1);
% $$$   patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
% $$$   jj = find(U > 1);
% $$$   patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
% $$$   drawnow;

    % visualization (somewhat slow)
    figure(1), clf,
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9]);
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    ii = find(U == 1);
    patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
    jj = find(U > 1);
    patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
    drawnow;

  % classify the DOFs

  adof = find(U);
  bdof_m = find(N*(U > 0) < neigh & U == 1);
  sdof = find(U > 1);
  sdof_m = find(N*(U > 1) < neigh & U > 1);
  idof = find(N*(U ~= 0) > 0 & U == 0);

  % DOFs which may proliferate
  pdof = find(U == 1);
  pdof = pdof(sqrt(P(1,pdof).^2+P(2,pdof).^2) <= rho);
  if isempty(pdof), pdof = zeros(0,1); end

  % "All DOFs" = adof + idof, like the "hull of adof"
  Adof = [adof; idof];
  % (after this adof is not used)
  
  % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
  % matrix. Determine also a local enumeration, eg. [1 2 3
  % ... numel(Adof)].
  Adof_ = (1:numel(Adof))';  
  [bdof_m_,sdof_,sdof_m_,idof_,pdof_] = ...
      map(Adof_,Adof,bdof_m,sdof,sdof_m,idof,pdof);

  % solve the currently active pressure Laplacian with direct injection
  % of Dirichlet BCs and sources
  if updLU
    La.X = L(Adof,Adof);
    Lai = fsparse(idof_,idof_,1,size(La.X));
    La.X = La.X-Lai*La.X+Lai;
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
    updLU = false;
  end
  Pr = full(fsparse(sdof_,1,1./dM(sdof),[size(La.X,1) 1]));
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr));

  % compute intensities of possible events

  % (1) moving DOFs

  % (1a) moving boundary DOFs
  grad = sum(N(bdof_m,idof),2).*Pr(bdof_m_);
  moveb = Drate*full(gradquotient*grad);

  % (1b) moving source DOFs
  [ii,jj_] = find(N(sdof_m,Adof));
  keep = find(U(Adof(jj_)) < 2);
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0),numel(sdof_m));
  moves = Drate*full(gradquotient*grad);
  
  % (2) cell proliferation
  prolif = lam*U(pdof)*(sum(U) < 1000);

  % waiting time until next event and the event itself
  intens = [moveb; moves; prolif];
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
  if tt+dt > Tend, dt = inf; end

  % record values as needed
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    IDsave(i+1:iend) = {ID};
    i = iend;
  end

  % reached the end (no final event)
  if isinf(dt)
    Usave(end) = {U};
    IDsave(end) = {ID};
    break;
  end

  if ix_ <= numel(moveb)
    % movement of a boundary (singly occupied) voxel
    ix = Adof(bdof_m_(ix_));
    n = find(N(ix,:));
    n = n(U(n) == 0);
    n = n(ceil(numel(n)*rand));

    % move from ix to n
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;

    % clear empty voxel
    ID(n,1) = ID(ix,1); ID(ix,1) = 0;
    updLU = true;
  elseif ix_ <= numel(moveb)+numel(moves)
    % movement of a cell in a doubly occupied voxel
    ix_ = ix_-numel(moveb);
    ix_ = sdof_m_(ix_);
    ix = Adof(ix_);

    jx_ = find(N(ix,Adof));
    jx_ = jx_(U(Adof(jx_)) < 2);
    rates = max(Pr(ix_)-Pr(jx_),0);
    m = find(cumsum(rates) > rand*sum(rates),1,'first');
    n = Adof(jx_(m));

    % execute event: move from ix to n and adjust format
    m = 1+(rand > 0.5); % 1 or 2
    if U(n) == 0
      % one cell moved into an empty voxel
      ID(n,1) = ID(ix,m);
      updLU = true;
    else
      % one cell moved into an already occupied voxel
      ID(n,2) = ID(ix,m);
    end
    if m == 1
      ID(ix,1) = ID(ix,2);
    end
    ID(ix,2) = 0;
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
  else
    % otherwise it was a cell division event
    ix_ = ix_-numel(moveb)-numel(moves);
    ix_ = pdof_(ix_);
    ix = Adof(ix_);

    U(ix) = U(ix)+1;
    % *** not implemented here: ID!
  end
  tt = tt+dt;
  report(tt,U,'');
end
report(tt,U,'done');

% display final frame
figure(1),
for i = numel(Usave)
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  axis([-1 1 -1 1]); axis square, axis off
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
  drawnow;
  title(sprintf('t = %f',tspan(i)));
end

return;

% postprocessing

% population appearance
M = struct('cdata',{},'colormap',{});
figure(1),
for i = 1:numel(Usave)
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
  axis([-1 1 -1 1]); axis square, axis off
  drawnow;
  title(sprintf('t = %f',tspan(i)));
  M(i) = getframe(gcf);
end

movie2gif(M,{M(1:2).cdata},'animations/growth.gif', ...
          'delaytime',0.1,'loopcount',0);
