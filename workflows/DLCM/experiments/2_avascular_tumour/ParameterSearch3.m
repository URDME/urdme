% Calibration of the Tumour/Proliferation/Deth/Degrading model.

% S. Engblom 2017-02-11
% D. Wilson 2017-09-13

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
Circ = zeros(Nvoxels);
for j = 1 : numel(xc)
  xp = xc(j);
  yp = yc(j);
  if ( xp < 1 || yp < 1 || xp > Nvoxels || yp > Nvoxels )
    continue
  end
  Circ(xp, yp, :) = 1;
end

% to calibrate: diffusive pressure rates
Drate1 = 0.01;     % into free matrix
%Drate2 = 25;      % into already visited matrix
Drate2 = 25;       % into already visited matrix (increased a bit)
%Drate3 = 0.01;    % into already occupied voxel
Drate3 = 0.01;     % into already occupied voxel (increased a bit)
Drate_ = [Drate1 Drate2; NaN Drate3];
% (hence Drate_(abs(U(j))+1,VU(j)+1) = Drate_(2*VU(j)+abs(U(j))+1) is
% the rate to move into voxel j)

% Drate1 <= Drate2 is the only thing that makes sense. And normalizing
% one of them to 1 sets the unit of time for the mechanics. Drate3 is
% less well defined as it partially is as result of the simplistic
% voxel-based description.

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% to calibrate: proliferation/death constants
cons = 0.0015;        % consumption of oxygen by cells
cutoff_prol = 0.65;  % (1 is the boundary source concentration)
r_prol = 0.125;         % rate of proliferation of singly occupied voxels
cutoff_die = 0.55;
r_die = 0.125;
r_degrade = 0.01;

K_birthmax = Inf;
birth_count=0;
% Cartesian mesh

h = 2/(Nvoxels-1);
xlayer = (-1:h:1)'; % [-1,1] x [-1,1]
ylayer = xlayer';
x = frepmat(xlayer,[1 numel(ylayer)]);
y = frepmat(ylayer,numel(xlayer));
P = [x(:)'; y(:)'];
extdof = find(Circ==1); % dofs for the sources at extreme boundary

% integration of pressure gradients over edges
D1 = h;
L1 = h;
gradquotient = L1/D1; % (= 1 for Cartesian mesh)

% triangulate the points
DT = delaunayTriangulation(P');
[V,R_] = voronoiDiagram(DT);
patchmax = max(cellfun('prodofsize',R_));
R = nan(size(R_,1),patchmax);
for i = 1:size(R,1), R(i,1:size(R_{i},2)) = R_{i}; end
% (allows for visualization by the patch-function)

% switch to PDE Toolbox mesh format [P,T]
P = DT.Points';
T = DT.ConnectivityList';
T = [T; ones(1,size(T,2))];

% assemble minus the Laplacian on this grid (ignoring BCs) as well as
% the Mass-matrix
[L,M] = assema(P,T,1,1,0);

% the (lumped) mass matrix gives the element volume
dM = full(sum(M,2));

% explicitly invert the lumped mass matrix and filter the diffusion matrix
[i,j,s] = find(L);
s = s./dM(i);
keep = find(s < 0);
i = i(keep); j = j(keep); s = s(keep);
L = fsparse(i,j,s,[Nvoxels Nvoxels].^2);
L = L+fsparse(1:Nvoxels^2,1:Nvoxels^2,-full(sum(L,2))');

% static neighbor matrix
N = fsparse(i,j,1,size(L));
neigh = full(sum(N,2));

% initial population, two cases:
init_case = 1;

if init_case == 1
  % (1) a small blob
  r = sqrt(P(1,:).^2+P(2,:).^2);
  ii = find(r < 0.2);
  U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
  
  % visit marker matrix: 1 for voxels who have been occupied
  VU = (U ~= 0);
elseif init_case == 2
  % (2) an approximation of the intended steady-state
  r = sqrt(P(1,:).^2+P(2,:).^2);
  ii = find(r < 0.35);
  U = fsparse(ii(:),1,1,[Nvoxels^2 1]);

  VU = (U ~= 0);

  % inner death radius about 10 cells
  ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.2);
  U0 = fsparse(ii(:),1,1,[Nvoxels^2 1]);

  % remove this
  %U = U-2*U0;
  % (unlike previous scripts these are now -1 and should degrade)
  U = U-U0;
else
  % case of U, UV already loaded
end

% visualization of this situation
% $$$ figure(1), clf,
% $$$ patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
% $$$       'EdgeColor','none');
% $$$ hold on,
% $$$ axis([-1 1 -1 1]); axis square, axis off
% $$$ ii = find(U == 1);
% $$$ patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
% $$$ jj = find(U == -1);
% $$$ patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[0 0 0]);
% $$$ title('Starting configuration: singles (green), dead (black)');
% $$$ drawnow;

% simulation interval
Tend = 500;
tspan = linspace(0,Tend,101);
report(tspan,'timeleft','init');

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;

tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);
while tt <= tspan(end)
  % classify the DOFs
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
  idof = find(Idof);
  
  Adof = [adof; idof];
  Adof_ = (1:numel(Adof))';
  [~,ix] = fsetop('ismember',bdof_m',Adof');
  bdof_m_ = Adof_(ix);
  [~,ix] = fsetop('ismember',sdof',Adof');
  sdof_ = Adof_(ix);
  [~,ix] = fsetop('ismember',sdof_m',Adof');
  sdof_m_ = Adof_(ix);
  [~,ix] = fsetop('ismember',idof1',Adof');
  idof1_ = Adof_(ix);
  [~,ix] = fsetop('ismember',idof2',Adof');
  idof2_ = Adof_(ix);
  [~,ix] = fsetop('ismember',idof',Adof');
  idof_ = Adof_(ix);
  [~,ix] = fsetop('ismember',adof',Adof');
  adof_ = Adof_(ix);
  
  if updLU
    % pressure Laplacian
    La.X = L(Adof,Adof);
    Lai = fsparse(idof_,idof_,1,size(La.X));
    La.X = La.X-Lai*La.X+Lai;
    [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

    % oxygen Laplacian
    OLa.X = L;
    OLai = fsparse(extdof,extdof,1,size(OLa.X));
    OLa.X = OLa.X-OLai*OLa.X+OLai;
    [OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

    updLU = false; % assume we can reuse
  end

  % RHS source term proportional to the over-occupancy and BCs
  Pr = full(fsparse(sdof_,1,1./dM(sdof), ...
                    [size(La.X,1) 1]));     % RHS first...
  Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

  % RHS source term proportional to the over-occupancy and BCs
  Oxy = full(fsparse([extdof; adof],1, ...
                     [ones(size(extdof)); ... 
                      -cons*full(max(U(adof),0)./dM(adof))], ...
                     [size(OLa.X,1) 1]));
  Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));

  % intensities of possible events

  % (1) moving boundary DOFs
  [ii,jj_] = find(N(bdof_m,Adof)); % neighbours...
  keep = find(U(Adof(jj_)) == 0);  % ...to move to
  ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
  % remove any possibly remaining negative rates
  grad = fsparse(ii,1,max(Pr(bdof_m_(ii))-Pr(jj_),0).* ...
                 Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                 numel(bdof_m));
  moveb = full(gradquotient*grad);

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
  birth = min(total_birth,K_birthmax)/total_birth * birth;
  birth(isnan(birth))=0; % as we get some 0/0 terms if total_birth=0;
 
  death = full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die));
  degrade = full(r_degrade*(U(Adof) == -1));
  
  intens = [moveb; moves; birth; death; degrade];
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

  % report back
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    i = iend;

    % monitor the maximum outlier cell:
    max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2))

    % the number of cells
    num_cells = sum(abs(U))

    % the rates
    inspect_rates = [sum(moveb) sum(moves) ...
                     sum(birth) sum(death) sum(degrade)]

    % visualization (slow things down a bit)
%     figure(2), clf,
%     patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%           'EdgeColor','none');
%     hold on,
%     axis([-1 1 -1 1]); axis square, axis off
%     ii = find(U == 1);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
%     ii = find(U == 2);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[1 0 0]);
%     ii = find(U == -1);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 0 0]);
%     title(sprintf('Time = %f, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
%     drawnow;
  end
  
  %%% visualisation real time
%   figure(2), clf,
%     patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%           'EdgeColor','none');
%     hold on,
%     axis([-1 1 -1 1]); axis square, axis off
%     ii = find(U == 1);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
%     ii = find(U == 2);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[1 0 0]);
%     ii = find(U == -1);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 0 0]);
%     title(sprintf('Time = %f, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
%     drawnow;
%     pause(0.0001)

  if ix_ <= numel(moveb)
    Ne.moveb = Ne.moveb+1;
    % movement of a boundary (singly occupied) voxel
    ix_ = bdof_m_(ix_);
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
  elseif ix_ <= numel(moveb)+numel(moves)
    Ne.moves = Ne.moves+1;
    % movement of a cell in a doubly occupied voxel
    ix_ = ix_-numel(moveb);
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
  elseif ix_ <= numel(moveb)+numel(moves)+numel(birth)
    Ne.birth = Ne.birth+1;
    % proliferation
    birth_count = birth_count+1;
    ix_ = ix_-numel(moveb)-numel(moves);
    ix = Adof(ix_);
    U(ix) = U(ix)+1;
  elseif ix_ <= numel(moveb)+numel(moves)+numel(birth)+numel(death)
    Ne.death = Ne.death+1;
    % death
    ix_ = ix_-numel(moveb)-numel(moves)-numel(birth);
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
    ix_ = ix_-numel(moveb)-numel(moves)-numel(birth)-numel(death);
    ix = Adof(ix_);
    U(ix) = 0;
    updLU = true; % boundary has changed
  end
  tt = tt+dt;
  report(tt,U,'');

  % update the visited sites
  VU = VU | U;
end

report(tt,U,'done');

% check that these are equal:
% $$$ sum(abs(Usave{1}))+Ne.birth-Ne.degrade
% $$$ full(sum(abs(Usave{end})))

return;

% population appearance
M = struct('cdata',{},'colormap',{});
figure(3), clf,
for i = 1:numel(Usave)
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
  hold on,
  axis([-1 1 -1 1]); axis square, axis off
  ii = find(Usave{i} == 1);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',graphics_color('bluish green'));
  ii = find(Usave{i} == 2);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',graphics_color('vermillion'));
  ii = find(Usave{i} == -1);
  patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 0 0]);
  title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
  drawnow;
  M(i) = getframe(gcf);
end

figure(4), clf
spsum  = @(U)(full(sum(abs(U))));
y = cellfun(spsum,Usave);
plot(tspan,y);
ylim([0 max(y)]);
xlabel('time')
ylabel('N cells')
