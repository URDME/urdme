% Delta-notch model in a growing population of cells.
%
%   This model runs the cell population physics from BASIC_TEST but
%   with a simple proliferation process and a delta-notch model at the
%   same time (both in continuous time). The purpose of the experiment
%   is to show how continuous-in-time processes may be easily coupled
%   together.

% S. Engblom 2017-12-28 (minor revision)
% S. Engblom 2017-08-30 (minor revision)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-27 (finalized physics)
% S. Engblom 2016-12-21 (updated the physics, changed to hex-cells)
% S. Engblom 2016-11-10

rng(12324); % repeatable results

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 60;

% diffusive pressure rate
Drate = 2; % (the simulation is rather sensitive to this parameter!)

% delta-notch timescale: 1 is "slow" compared to mechanics, 50 or
% above is "fast"
Tdn = 1;

% fetch discretization
[P,E,T,gradquotient] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);
neigh = 6;

% birth process
lam = 1;    % (initial) proliferation rate
rho = 0.10; % radius of proliferation

% initial population
ii = find(sqrt(sum(P.^2)) <= 2*rho);
U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
% sparse solution vector U(i) = {0,1,2} cells for voxel i

% delta-notch DOFs
Delta = fsparse(ii(:),1,rand(numel(ii),1),[Nvoxels^2 2]);
Notch = fsparse(ii(:),1,rand(numel(ii),1),[Nvoxels^2 2]);
% Format: these are a sparse solutions vector over the full mesh
% [P,E,T], where the second column is unused until cells proliferate,
% a filled row then means the corresponding voxel contains two cells
% (which is the maximum).

% delta-notch model
a = 0.01;
b = 100;
v = 1;
k = 2;
h = 2;

F = @(x)(x.^k./(a+x.^k));
G = @(x)(1./(1+b*x.^h));

sol = 1; % Wiener SDE interpretation on/off
Vol = 16; iVol = 1/Vol; % (imagined) system volume

% max time-step in forward Euler
maxdt = 0.5/Tdn;

% simulation interval
Tend = 100;
% solution recorded at this grid:
tspan = linspace(0,Tend,61);

report(tspan,U,'init');

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Dsave = cell(1,numel(tspan));
Nsave = cell(1,numel(tspan));

tt = tspan(1);
Usave{1} = U;
Dsave{1} = Delta;
Nsave{1} = Notch;
i = 1;
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
  
  % local enumeration
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

  % use the waiting time dt to update Delta-Notch (forward Euler)
  ddt = min(Tend-tt,dt);
  Nsteps = max(ceil(ddt/maxdt),1);
  ddt = ddt/Nsteps;
  for j = 1:Nsteps
    % record values as needed
    if tspan(i+1) < tt+ddt
      iend = i+find(tspan(i+1:end) < tt+ddt,1,'last');
      Usave(i+1:iend) = {U};
      Dsave(i+1:iend) = {Delta};
      Nsave(i+1:iend) = {Notch};
      i = iend;
    end

    % average incoming Delta from neighbours
    Delta_ = sum(N*Delta,2);
    M = N*U;
    M(sdof) = M(sdof)+1; % (one extra in the same voxel)
 
    if sol == 1
      % Forward Euler step for the delta-notch model
      %  D' = v (G(N)-D)
      %  N' = F(D_mean)-N
      Delta(adof,1) = Delta(adof,1)+ddt*Tdn*v*(G(Notch(adof,1))-Delta(adof,1));
      Notch(adof,1) = Notch(adof,1)+ ...
          ddt*Tdn*(F((Delta_(adof)+Delta(adof,2))./max(M(adof),1))-Notch(adof,1));
      % doubly occupied cells:
      Delta(sdof,2) = Delta(sdof,2)+ddt*Tdn*v*(G(Notch(sdof,2))-Delta(sdof,2));
      Notch(sdof,2) = Notch(sdof,2)+ ...
          ddt*Tdn*(F((Delta_(sdof)+Delta(sdof,1))./max(M(sdof),1))-Notch(sdof,1));
    else
      % The corresponding Wiener SDE becomes
      %  dD = v (G(N)-D) dt + v^(1/2) [G(N)^(1/2) dW_1 - D^(1/2) dW_2]
      %  dN = [F(D_mean)-N] dt + F(D_mean)^(1/2) dW_3 - N^(1/2) dW_4
      error('Not implemented.');
    end

    tt = tt+ddt;
  end
  % reached the end (no final event)
  if isinf(dt)
     Usave(end) = {U};
     Dsave(end) = {Delta};
     Nsave(end) = {Notch};
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
    Delta(n,1) = Delta(ix,1); Delta(ix,1) = 0;
    Notch(n,1) = Notch(ix,1); Notch(ix,1) = 0;
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
      Delta(n,1) = Delta(ix,m);
      Notch(n,1) = Notch(ix,m);
      updLU = true;
    else
      % one cell moved into an already occupied voxel
      Delta(n,2) = Delta(ix,m);
      Notch(n,2) = Notch(ix,m);
    end
    if m == 1
      Delta(ix,1) = Delta(ix,2);
      Notch(ix,1) = Notch(ix,2); 
    end
    Delta(ix,2) = 0;
    Notch(ix,2) = 0;
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
  else
    % otherwise it was a cell division event
    ix_ = ix_-numel(moveb)-numel(moves);
    ix_ = pdof_(ix_);
    ix = Adof(ix_);

    U(ix) = U(ix)+1;
    % share Delta-Notch between the two cells
    r = rand; Delta(ix,2) = r*Delta(ix,1); Delta(ix,1) = (1-r)*Delta(ix,1);
    r = rand; Notch(ix,2) = r*Notch(ix,1); Notch(ix,1) = (1-r)*Notch(ix,1);
  end
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
  patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'));
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'));
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
  patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'));
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'));
  axis([-1 1 -1 1]); axis square, axis off
  drawnow;
  title(sprintf('t = %f',tspan(i)));
  M(i) = getframe(gcf);
end

movie2gif(M,{M(1:2).cdata},'animations/delta_notch_population.gif', ...
          'delaytime',0.1,'loopcount',0);

% delta-notch appearance
M2 = struct('cdata',{},'colormap',{});
figure(2),
for i = 1:numel(Usave)
  clf,
  patch('Faces',R,'Vertices',V, ...
        'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  patch('Faces',R(Nsave{i}(:,1) > 0.5 & Usave{i} > 0,:), ...
        'Vertices',V,'FaceColor',[0 0 0]);
  patch('Faces',R(Nsave{i}(:,1) <= 0.5 & Usave{i} > 0,:), ...
        'Vertices',V,'FaceColor',[0.6 0.6 0.6]);
  axis([-1 1 -1 1]); axis square, axis off
  drawnow;
  M2(i) = getframe(gcf);
end

if Tdn == 1
  movie2gif(M2,{M2(1:2).cdata},'animations/delta_notch.gif', ...
            'delaytime',0.1,'loopcount',0);
else
  movie2gif(M2,{M2(1:2).cdata},'animations/delta_notch_fast.gif', ...
            'delaytime',0.1,'loopcount',0);
end

% a few such frames
j = 2;
for i = ceil([0.05 0.4 0.7 1]*numel(Usave))
  j = j+1;
  figure(j), clf,
  patch('Faces',R,'Vertices',V, ...
        'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  patch('Faces',R(Nsave{i}(:,1) > 0.5 & Usave{i} > 0,:), ...
        'Vertices',V,'FaceColor',[0 0 0]);
  patch('Faces',R(Nsave{i}(:,1) <= 0.5 & Usave{i} > 0,:), ...
        'Vertices',V,'FaceColor',[0.6 0.6 0.6]);
  axis([-1 1 -1 1]); axis square, axis off
  if j == 3
    z = exp(1i*2*pi*linspace(0,1,40))*rho;
    plot(z,'w','linewidth',2);
  end
  drawnow;
  
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 380 260]);
  drawnow;
% $$$   % uncomment to save:
% $$$   set(gcf,'inverthardcopy','off')
% $$$   if Tdn == 1
% $$$     print('-depsc',sprintf('figures/delta_notch%d.eps',j-2));
% $$$   else
% $$$     print('-depsc',sprintf('figures/delta_notch%d_fast.eps',j-2));
% $$$   end
end
