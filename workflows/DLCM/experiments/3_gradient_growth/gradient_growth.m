% Simulation of a population of cells responding to a chemical gradient.

% S. Engblom 2017-12-28 (minor revision)
% S. Engblom 2016-11-09 (revision)
% S. Engblom 2016-07-05 (minor revision)
% S. Engblom 2016-05-01

%% discretization (used ndop rather than dt_operators)

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 100;

% assemble minus the Laplacian on this grid (ignoring BCs)
L = -ndop([0 1 0; 1 -4 1; 0 1 0],[2 2],[Nvoxels Nvoxels]);
% gradient matrices in the x/y-directions
X = ndop([0 0 0; -1 0 1; 0 0 0],[2,2],[Nvoxels Nvoxels]);
Y = ndop([0 -1 0; 0 0 0; 0 1 0],[2,2],[Nvoxels Nvoxels]);

%% gradient of the SLIT chemical

% parameters for the diffusion and degradation of the chemical
Q = 2; k = 0.1; Ds = 50; lambda = sqrt(k/Ds);
xpos = 1:Nvoxels;
sx = Q/(2*Ds*lambda)*exp(-lambda*(Nvoxels-xpos));
S = repmat(sx,Nvoxels,1);
Sx = spdiags(X*S(:),0,Nvoxels^2,Nvoxels^2);
Sl = spdiags(L*S(:),0,Nvoxels^2,Nvoxels^2);

%% advection-diffusion operator
chi = -1250;
A = L+chi*Sl+0.25*chi*Sx*X; % we ignore Sy as SLIT is constant in y

% neighbor matrix (used to select the voxel [left above below right]
N = ndop([0 1 0; 1 0 1; 0 1 0],[2 2],[Nvoxels Nvoxels]);

% initial populations of cells
XX = repmat(1:Nvoxels,Nvoxels,1);
YY = XX';
U = 3*((XX-Nvoxels/2).^2+(YY-Nvoxels/2).^2 < 100);

% simulation interval
Tend = 1000;
% solution recorded at this grid:
tspan = linspace(0,Tend,100);

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
while tt <= tspan(end)
  % classify the Degrees-Of-Freedoms (DOFs for short)
  adof = find(U); % active DOFs: occupied voxels
  bdof = find(N*(U(:) > 0) < 4 & U(:)); % boundary DOFs

  % The above will be enumerated within U, an imagined
  % Nvoxels-by-Nvoxels sparse matrix. Determine also a local
  % enumeration, eg. [1 2 3 ... numel(adof)].
  adof_ = (1:numel(adof))';
  [~,ix] = fsetop('ismember',bdof',adof');
  bdof_ = adof_(ix);

  % the currently active Laplacian
  La = A(adof,adof);
  % (note: this indexing/truncation implies a homogeneous Dirichlet
  % boundary conditions in the layer of voxels just outside the
  % boundary)

  % Solve the pressure equation with RHS source terms proportional to
  % the over-occupancy.
  P = La\(full(U(adof))-1);
  % (note: the -1 is to offset the 0 pressure of the free matrix)

  % Movement rates are proportional to minus the gradient (also with
  % homogeneous Dirichlet boundary conditions in the layer of voxels
  % just outside the boundary).

  % only boundary DOFs can move
  Press = zeros(Nvoxels);
  Press(adof) = P;
  ii = find(U);
  intens = sqrt((Press(ii(bdof_)+1) - Press(ii(bdof_)-1)).^2 + ...
                (Press(ii(bdof_)+Nvoxels) - Press(ii(bdof_)-Nvoxels)).^2);

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
  % (now ix_ points to the boundary DOF which will move)

  % record values as needed
  if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    i = iend;
  end

  % find DOF of event in global enumeration
  ix = adof(bdof_(ix_));

  % neighbors of boundary DOF ix... which are empty... which is the one
  % selected at random
  n = find(N(ix,:));
  n = n(U(n) == 0);
  n = n(ceil(numel(n)*rand));

  % execute event: move from ix to n
  U(ix) = U(ix)-1;
  U(n) = U(n)+1;

  tt = tt+dt;
  report(tt,U,'');
end
report(tt,U,'done');

return;

% do the above, e.g., 1000 times to compute an average:
Store = zeros(Nvoxels);
for i = 1:1000
  gradient_growth;
  Store = Store + U;
end
Store = Store/1000;
% uncomment to save to file:
% save experiments/3_gradient_growth/1000_realisations

return;

% graphics for this
load 1000_realisations

% Cartesian mesh
[P,E,T] = basic_mesh(1,Nvoxels);
P = P([2 1],:);
[V,R] = mesh2dual(P,E,T,'voronoi');

h = patch('Faces',R,'Vertices',V,'FaceVertexCData',Store(:), ...
          'FaceColor','flat','EdgeColor','none');
hold on,
axis([-sqrt(3)/2 sqrt(3)/2 -1 1]);
axis equal

z = (sqrt(10)/15)*exp(2i*pi*linspace(0,1,40));
plot(z,'k','LineWidth',1);

set(gca,'xtick',[],'ytick',[],'Visible','off');
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
colorbar

% uncomment to print to file:
%print -depsc gradient_growth.eps
