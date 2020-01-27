% Calibration of the Avascular Tumour Model .

% S. Engblom 2017-02-11
% D. Wilson 2017-09-13

% Order of search:
%
% (1) decide where the steady-state should be
%
% ==> Decision: a tumour of radius 25, with (inner-) death radius 10.
%
% (2) solve the Oxygen equation on this steady-state and decide
% suitable values of (cons,cutoff_prol,cutoff_die) such that (i)
% proliferation only occurs in the (say) two outer layers of cells,
% and (ii), that death occurs in the (say) two layers just outside the
% death radius
%
% ==> Calibration: (cons,cutoff_prol,cutoff_die) = (0.01,0.9,0.05).
% This is the outcome of ParameterSearch1 (the current script).
%
% (3) next, pretend that r_prol ~ r_die = <something whatever> and
% that r_degrade = inf such that dead cells degrade immediately
%
% (4) this means that the inner death radius is actually a hole, so we
% can now solve the pressure equation with this hole in place and
% sources in the two outer layers of the tumour to decide suitable
% values of Drate such that the sum of the pressure gradients for all
% cells to go into the hole is (say) 100 times bigger than the sum of
% the pressure gradients to go into the free matrix
%
% ==> Calibration: Drate1 = 0.01, Drate2 = 25, Drate3 = <not decided>.
% This is the outcome of the script ParameterSearch2.
%
% (5) from here we go to the script ParameterSearch3 and fiddle a bit
% with (r_prol,r_die,r_degrade) + initial population and possibly the
% other constants to get a quasi-static behavior.

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
Circ=zeros(Nvoxels);
for j = 1 : numel(xc)
  xp = xc(j);
  yp = yc(j);
  if ( xp < 1 || yp < 1 || xp > Nvoxels || yp > Nvoxels )
    continue
  end
  Circ(xp, yp, :) = 1;
end

% to calibrate: diffusive pressure rates
Drate1 = 0.01;  % into free matrix
Drate2 = 100;   % into already visited matrix
Drate3 = 0.01;  % into already occupied voxel
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
cons = 0.0125;        % consumption of oxygen by cells
cutoff_prol = 0.65;  % (1 is the boundary source concentration)
r_prol = 0.1;       % the rate of proliferation of voxels where U==1
cutoff_die = 0.55;
r_die = 0.01;
r_degrade = 0.0;

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
ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.35);
U = fsparse(ii(:),1,1,[Nvoxels^2 1]);

% outer radius about 25 cells
radius = sqrt(full(sum(U))/pi)

% visit marker matrix: 1 for voxels who have been occupied
VU = (U > 0);

% inner death radius about 10 cells
ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.2);
U0 = fsparse(ii(:),1,1,[Nvoxels^2 1]);
radius0 = sqrt(full(sum(U0))/pi)

% remove this
U = U-U0;

% classify the DOFs
adof = find(U); % all filled voxels
bdof_m = find(N*(U ~= 0) < neigh & abs(U) == 1); % the singularly occupied voxels on the boundary
sdof = find(U > 1); % the voxels with 2 cells in them
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
Oxy = full(fsparse([extdof; adof],1, ...
                   [ones(length(extdof),1); ... 
                    -cons*full(max(U(adof),0)./dM(adof))], ...
                   [size(OLa.X,1) 1]));
Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));

% visualize the oxygen concentration
figure(1), clf
plot3(P(1,Adof),P(2,Adof),Oxy(Adof),'b.');
% (note: the max(.,0) is the concentration the cells "see")
hold on,

% time to calibrate these:
%cutoff_prol = 0.9;
%cutoff_die = 0.05;
ii = find(Oxy(Adof) > cutoff_prol);
jj = find(Oxy(Adof) < cutoff_die);
plot3(P(1,Adof(ii)),P(2,Adof(ii)),max(Oxy(Adof(ii)),0),'ro');
plot3(P(1,Adof(jj)),P(2,Adof(jj)),max(Oxy(Adof(jj)),0),'mo');
legend('cells','proliferating','dying');
title('Oxygen');
