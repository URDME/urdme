%GROWTH
%
%   Growth of domain when all cells under a high enough oxygen
%   concentration are allowed to proliferate. Exports a full event
%   history (Event,Etime,Etype).

% S. Engblom 2018-01-25 (new context)
% S. Engblom 2017-08-30 (minor revision)
% S. Engblom 2016-12-28 (reuse of factorization)
% S. Engblom 2016-12-27 (new physics)
% S. Engblom 2016-11-30

% diffusive pressure rate
Drate = 1;

% initial population
Nvoxels = size(U,1);
if ~exist('Usave','var')
  [~,ii] = min(sqrt(P(1,:).^2+P(2,:).^2));
  U = fsparse(ii(:),1,1,[Nvoxels 1]);
else
  U = Usave{end};
  warning('Restart from previous values.');
end

% oxygen consumption rate per cell
gamma = 0.05;

% proliferation: rate and cutoff
lam = 0.5;
prolif_cutoff = 0.5; % (1 is the boundary source concentration)

% simulation interval
if ~exist('Tend','var')
  warning('Non-existing variable (Tend). Default used.');
  Tend = 40;
end

% solution recorded at this grid:
tspan = linspace(0,Tend,51);
report(tspan,'timeleft','init');

% event recorder in expanding arrays
Etime = zeros(1,0);
Eii = zeros(1,0);
Ejj = zeros(1,0); ejj = 0; % counter
Evv = zeros(1,0);
Etype = zeros(1,0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;

tt = tspan(1);
i = 1;
neigh = 6;
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
  
  % classify the DOFs
  adof = find(U);
  bdof_m = find(N*(U > 0) < neigh & U == 1);
  sdof = find(U > 1);
  sdof_m = find(N*(U > 1) < neigh & U > 1);
  udof = find(U == 1); % may proliferate
  idof = find(N*(U ~= 0) > 0 & U == 0);
  Adof = [adof; idof];

  % local enumeration
  Adof_ = (1:numel(Adof))';
  Adof_ = (1:numel(Adof))';  
  [adof_,bdof_m_,sdof_,sdof_m_,udof_,idof_] = ...
      map(Adof_,Adof,adof,bdof_m,sdof,sdof_m,udof,idof);

  if updLU
    % pressure/oxygen Laplacian - reused factorization
    La.X = L(Adof,Adof);
    Lai = fsparse(idof_,idof_,1,size(La.X));
    La.X = La.X-Lai*La.X+Lai;
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
    updLU = false; % assume we can reuse
  end

  % pressure sources: over-occupied voxels
  Pr = full(fsparse(sdof_,1,1./dM(sdof),[size(La.X,1) 1]));
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr));
  % oxygen sources: boundary voxels, sinks: cells
  Oxy = full(fsparse([idof_; adof_],1, ...
                     [ones(numel(idof_),1); -gamma*full(U(adof)./dM(adof))], ...
                     [size(La.X,1) 1]));
  Oxy(La.q) = La.U\(La.L\(La.R(:,La.p)\Oxy));

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

  % (2) proliferation
  if sum(U) >= 1000
    prolif = zeros(0,1);
  else
    pdof_ = udof_(Oxy(udof_) > prolif_cutoff);
    prolif = lam*Oxy(pdof_);
  end

  % waiting time until next event and the event itself
  intens = [moveb; moves; prolif];
  a0 = sum(intens);
  dt = -reallog(rand)/a0;
  rnd = rand*a0;
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
    if iend == numel(tspan)
      break;
    end
  end

  Etime = [Etime tt+dt];
  ejj = ejj+1;
  if ix_ <= numel(moveb)
    % movement of a boundary (singly occupied) voxel
    ix = Adof(bdof_m_(ix_));
    n = find(N(ix,:));
    n = n(U(n) == 0);
    n = n(ceil(numel(n)*rand));

    % execute event: move from ix to n
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
    Eii = [Eii ix n];
    Ejj = [Ejj ejj ejj];
    Evv = [Evv -1 1];
    Etype = [Etype 1];
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

    % execute event: move from ix to n
    if U(n) == 0
      updLU = true;
    end
    U(ix) = U(ix)-1;
    U(n) = U(n)+1;
    Eii = [Eii ix n];
    Ejj = [Ejj ejj ejj];
    Evv = [Evv -1 1];
    Etype = [Etype 2];
  else
    % proliferation event
    ix_ = ix_-numel(moveb)-numel(moves);
    ix_ = pdof_(ix_);
    ix = Adof(ix_);
    U(ix) = 2;
    Eii = [Eii ix];
    Ejj = [Ejj ejj];
    Evv = [Evv 1];
    Etype = [Etype 3];
  end
  tt = tt+dt;
  report(tt,U,'');
end

% (very) sparse event matrix
Event = fsparse(Eii,Ejj,Evv,[size(U,1) size(Etime,2)]);

% check:
% norm(Usave{1}+sum(Event,2)-Usave{end},1) == 0

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

% population appearance
figure(1),
for i = 1:numel(Usave)
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
  hold on,
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0],'EdgeColor','none');
  jj = find(Usave{i} > 1);
  patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0],'EdgeColor','none');
  axis([-1 1 -1 1]); axis square, axis off
  drawnow;
  title(sprintf('t = %f',tspan(i)));
end
