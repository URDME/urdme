%SCHEMATICS Draws schematics of cell mechanics.
%
%   Schematics of PDE-based cell mechanics, draws an explanatory
%   schematic picture.

% S. Engblom 2018-01-25

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

neigh = full(max(sum(N,2)));

% population
ii = find(sqrt(sum(P.^2)) <= 0.1);
ii = [ii find(sqrt(sum(tsum(P,-[-0.2 -0.2],[1 2],[3 1]).^2)) <= 0.05)];
jj = find(sqrt(sum(P.^2)) <= 0.3);

% assemble this
U = fsparse([ii(:); jj(:)],1,1,[Nvoxels^2 1]);

figure(4), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,
axis([-0.4 0.2 -0.4 0.2]);
axis square, axis off
ii = find(U == 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'));
jj = find(U > 1);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'));

% "The Method"
adof = find(U);
bdof_m = find(N*(U > 0) < neigh & U == 1);
sdof = find(U > 1);
sdof_m = find(N*(U > 1) < neigh & U > 1);
idof = find(N*(U ~= 0) > 0 & U == 0);
Adof = [adof; idof];

Adof_ = (1:numel(Adof))';  
[bdof_m_,sdof_,sdof_m_,idof_] = ...
    map(Adof_,Adof,bdof_m,sdof,sdof_m,idof);

La = L(Adof,Adof);
Lai = fsparse(idof_,idof_,1,size(La));
La = La-Lai*La+Lai;

Pr = full(La\fsparse(sdof_,1,1./dM(sdof),[size(La,1) 1]));

grad = sum(N(bdof_m,idof),2).*Pr(bdof_m_);
moveb = Drate*full(gradquotient*grad);

% visualize this
x = zeros(1,0); y = zeros(1,0);
u = zeros(1,0); v = zeros(1,0);
k_ = 0;
for i = bdof_m'
  k_ = k_+1;
  for j = idof(find(N(i,idof)))'
    x = [x P(1,i)];
    y = [y P(2,i)];
    u = [u (P(1,j)-P(1,i))*Drate*gradquotient*Pr(bdof_m_(k_))];
    v = [v (P(2,j)-P(2,i))*Drate*gradquotient*Pr(bdof_m_(k_))];
  end
end

[i,j_] = find(N(sdof_m,Adof));
keep = find(U(Adof(j_)) < 2);
i = i(keep); j_ = j_(keep);
grad = fsparse(i,1,max(Pr(sdof_m_(i))-Pr(j_),0),numel(sdof_m));
moves = Drate*full(gradquotient*grad);

% visualize this
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
      u = [u (P(1,Adof(j_))-P(1,i))*Drate*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
    end
  end
end
quiver(x,y,u,v,'k');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
drawnow;
% uncomment to save:
% $$$ print -depsc ../../figures/basic_schematics.eps
