% Mechanics explanation.
%
%   This script 'explains' the cell mechanics by creating a sample
%   population of cells and then determining all movement rates and
%   visualizing them. The use of this could be for quick prototyping
%   or pedagogical purposes.

% S. Engblom 2017-12-26 (revision)
% S. Engblom 2017-08-29 (minor revision)
% S. Engblom 2017-01-07

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 40;

% diffusive pressure rates
Drate1 = 0.1;    % into free matrix
Drate2 = 1;      % into already visited matrix
Drate3 = 0.1; % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];
% (hence Drate_(U(j)+1,VU(j)+1) = Drate_(2*VU(j)+U(j)+1) is the rate
% to move into voxel j). This works due to 2*VU(j)+U(j)+1 == {1 --> Drate1, 
% 2 --> NaN, 3 --> Drate2, 4 --> Drate3}

% Drate1 <= Drate2 is the only thing that makes sense. And normalizing
% one of them to 1 sets the unit of time for the mechanics. Drate3 is
% less well defined as it partially is a result of the simplistic
% voxel-based description.

% boundary conditions: free matrix BC1, already visited matrix BC2 (0:
% no pressure, +: must overcome a pressure, -: cells get "sucked" out)
BC1 = 0;
BC2 = 0;

% BC1 >= BC2 is the only thing that makes sense. And normalizing one
% of them to 0 sets the pressure unit.

% fetch discretization (mesh_type = 1 or 2)
if ~exist('mesh_type','var'), error('Must define mesh_type.'); end
[P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);

% initial population
ii = find(abs(P(2,:)+1) < 1);
U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
% visit marker matrix: 1 for voxels who have been occupied
VU = (U > 0);

ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.5);
U(ii) = 0; % a whole

ii = find(sqrt((P(1,:)-0.7).^2+((P(2,:)+0.3).^2)) < 0.1);
ii = [ii find(sqrt((P(1,:)+0.75).^2+((P(2,:)+0.15).^2)) < 0.07)];
U(ii) = 2; % doubly occupied blobs

ii = find(abs(P(2,:)+1) < 0.35);
U(ii) = 0; % an empty bottom strip

% classify the DOFs

neigh = full(sum(N,2));
adof = find(U);
bdof_m = find(N*(U > 0) < neigh & U == 1);
sdof = find(U > 1);
sdof_m = find(N*(U > 1) < neigh & U > 1);
Idof = (N*(U ~= 0) > 0 & U == 0);
idof1 = find(Idof & ~VU); % "external" BC1
idof2 = find(Idof & VU);  % "internal" BC2
idof = find(Idof);

Adof = [adof; idof];
Adof_ = (1:numel(Adof))';  
[bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_] = ...
    map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof);

% pressure Laplacian - factorization is reused
La.X = L(Adof,Adof);
Lai = fsparse(idof_,idof_,1,size(La.X));
La.X = La.X-Lai*La.X+Lai;
[La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');

% RHS source term proportional to the over-occupancy and BCs
Pr = full(fsparse([sdof_; idof1_; idof2_],1, ...
                  [1./dM(sdof); ...
                   BC1*ones(size(idof1_)); BC2*ones(size(idof2_))], ...
                  [size(La.X,1) 1]));    % RHS first...
Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

% intensities of possible events

% (1) moving boundary DOFs
% this code assumes homogeneous Dirichlet BCs + a single Drate:
% $$$   grad = sum(N(bdof_m,idof),2).*Pr(bdof_m_);
% $$$   moveb = Drate*full(gradquotient*grad);
% this code takes care of arbitrary Dirichlet BCs + different Drates:
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
               Drate_(2*VU(Adof(jj_))+U(Adof(jj_))+1), ...
               numel(sdof_m));
moves = full(gradquotient*grad);

 % pressure plot
 figure(3), clf,
 Pr_ = full(U); Pr_(Adof) = Pr;
 trisurf(T(1:3,:)',P(1,:),P(2,:),Pr_);

% visualize rates
figure(mesh_type), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,
axis equal, axis off
ii = find(U == 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'));
jj = find(U > 1);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'));
title(sprintf('BC = [%.0f,%.0f], Drate = [%.1e,%.1e,%.1e]', ...
                BC1,BC2,Drate1,Drate2,Drate3));

% visited boundary voxel or not
plot(P(1,idof1),P(2,idof1),'ko');
plot(P(1,idof2),P(2,idof2),'k.');

% compute all "moveb"-rates
x = zeros(1,0); y = zeros(1,0);
u = zeros(1,0); v = zeros(1,0);
k_ = 0;
for i = bdof_m'
  k_ = k_+1;
  jj_ = find(N(i,Adof));
  keep = find(U(Adof(jj_)) < 1);
  jj_ = jj_(keep);
  for j_ = jj_
    if Pr(bdof_m_(k_)) > Pr(j_)
      x = [x P(1,i)];
      y = [y P(2,i)];
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient*(Pr(bdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient*(Pr(bdof_m_(k_))-Pr(j_))];
    end
  end
end

% compute all "moves"-rates
k_ = 0;
for i = sdof_m'
  k_ = k_+1;
  jj_ = find(N(i,Adof));
  keep = find(U(Adof(jj_)) < 2);
  jj_ = jj_(keep);
  for j_ = jj_
    if Pr(sdof_m_(k_)) > Pr(j_)
      x = [x P(1,i)];
      y = [y P(2,i)];
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
    end
  end
end
quiver(x,y,u,v);
drawnow;
