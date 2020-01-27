% Calibration of the Tumour/Proliferation/Deth/Degrading model.

% S. Engblom 2017-02-11
% D. Wilson 2017-09-13

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 101; % odd so the BC for oxygen can by centered

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
Drate1 = 0.01;    % into free matrix
Drate2 = 10;      % into already visited matrix
Drate3 = 0.01;    % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];
% (hence Drate_(abs(U(j))+1,VU(j)+1) = Drate_(2*VU(j)+abs(U(j))+1) is
% the rate to move into voxel j)

% Drate1 <= Drate2 is the only thing that makes sense. And normalizing
% one of them to 1 sets the unit of time for the mechanics. Drate3 is
% less well defined as it partially is as result of the simplistic
% voxel-based description.

% boundary conditions
OBC1 = 1; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% to calibrate: proliferation/death constants
cons = 0.0025;        % consumption of oxygen by cells
cutoff_prol = 0.725;  % (1 is the boundary source concentration)
r_prol = 0.1;       % the rate of proliferation of voxels where U==1
cutoff_die = 0.65;
r_die = 0.01;
r_degrade = 0.0;

K_birthmax = 10;

% Cartesian mesh

h = 2/(Nvoxels-1);
xlayer = (-1:h:1)'; % [-1,1] x [-1,1]
ylayer = xlayer';
x = frepmat(xlayer,[1 numel(ylayer)]);
y = frepmat(ylayer,numel(xlayer));
P = [x(:)'; y(:)'];

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

% initial population
r = sqrt(P(1,:).^2+P(2,:).^2);
ii = find(r < 0.35);
U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
% extra layer of about two cells which are doubly occupied:
ii = find(0.32 < r & r < 0.35);
U = U+fsparse(ii(:),1,1,[Nvoxels^2 1]);

% visit marker matrix: 1 for voxels who have been occupied
VU = (U > 0);

% inner death radius about 10 cells
ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.2);
U0 = fsparse(ii(:),1,1,[Nvoxels^2 1]);

% remove this
U = U-U0;

% visualization of this situation
figure(1), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
      'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
ii = find(U == 1);
patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);
jj = find(U > 1);
patch('Faces',R(jj,:),'Vertices',V,'FaceColor',[1 0 0]);
title('Test population: sources (red), singles (green)');
drawnow;

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
extdof = find(Circ==1);
  
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

% RHS source term proportional to the over-occupancy and BCs
Pr = full(fsparse(sdof_,1,1./dM(sdof), ...
                  [size(La.X,1) 1]));     % RHS first...
Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

% RHS source term proportional to the over-occupancy and BCs
Oxy = full(fsparse([extdof; adof_],1, ...
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

% visualization of this situation
figure(2), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
      'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
ii = find(U);
patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0 1 0]);

patch('Faces',R(bdof_m,:),'Vertices',V,'FaceColor',[0 0 1]);
patch('Faces',R(sdof_m,:),'Vertices',V,'FaceColor',[1 0 0]);
title('Moving sources (red), moving singles (blue)');
drawnow;

% should be about 100 times:
moveb_to_moves = sum(moveb)/sum(moves)
