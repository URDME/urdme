% Modeling explanation.
%
%   This script draws an illustrative schematic picture explaining the
%   numerical model.

% E. Blom 2023-22-09 (minor revision)
% S. Engblom 2017-12-26 (revision)
% S. Engblom 2017-08-29 (minor revision)
% S. Engblom 2017-01-30

load experiments/petri_mesh/petri_geom % load [G,B]

% create mesh
[P,E,T] = initmesh(G,'hmax',0.35);   % uniform mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);

% each node in the mesh is a (midpoint of a) voxel
Nvoxels = size(P,2);

% radius of each voxel
Radius = R;
Radius(isnan(Radius)) = 1; % map all NaN's to the infinity vertex
Radius = sqrt(min(sum(tsum(reshape(V(Radius,:),Nvoxels,size(R,2),2),-P, ...
                           [1 3 2],[2 1]).^2,2),[],3));

% population
% $$$ ii = find(max(abs(P),[],1) <= 0.7);
% $$$ U = fsparse(ii(:),1,1,[Nvoxels 1]);
% $$$ jj = find(max(abs(P),[],1) <= 0.1);
% $$$ U = U+fsparse(jj(:),1,1,[Nvoxels 1]);
ii = find(max(abs(P),[],1) <= 0.6);
U = fsparse(ii(:),1,1,[Nvoxels 1]);
jj = find(max(abs(P),[],1) <= 0.1);
U = U+fsparse(jj(:),1,1,[Nvoxels 1]);

%% Cell representation in voxels
figure(1), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off

ii = find(U == 1);
theta = linspace(0,2*pi);
for i = ii'
  z = complex(P(1,i),P(2,i))+Radius(i)*exp(1i*theta);
  plot(z,'k:','LineWidth',1.5);
  z = complex(P(1,i),P(2,i))+0.2*Radius(i)*exp(1i*theta);
  plot(z,'k:','LineWidth',1.5);
end

jj = find(U > 1);
for j = jj'
  % two ellipses make up the crammed up cells
  for displ = [-0.05 0.05; 0.05 -0.03; 0.85 0.9]
    for scale = [0.8 0.15]
      x1 = displ(1)+P(1,j)-scale*Radius(j);
      y1 = displ(2)+P(2,j)-scale*Radius(j);
      x2 = displ(1)+P(1,j)+scale*Radius(j);
      y2 = displ(2)+P(2,j)+scale*Radius(j);
      e = displ(3);
      a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
      b = a*sqrt(1-e^2);
      X = a*cos(theta);
      Y = b*sin(theta);
      w = atan2(y2-y1,x2-x1);
      x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
      y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
      plot(x,y,'k','LineWidth',1.5);
    end
  end
end

jj = find(U > 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'),'FaceAlpha',0.6);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'),'FaceAlpha',0.6);
axis([-0.25 1 -0.25 1]);

drawnow;

%% a second time, but with "connecting springs"
figure(2), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off

ii = find(U == 1);
jj = find(U > 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'),'FaceAlpha',0.6);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'),'FaceAlpha',0.6);
axis([-0.25 1 -0.25 1]);

ii = find(U);
nrounds = 6;
theta = linspace(0,2*pi*nrounds,100*nrounds);
mask = ones(size(theta));
mask([1:200 end-199:end]) = [linspace(0,1,200) 1-linspace(0,1,200)];
mask = mask.^6;
for i = ii'
  jj = find(N(i,:) & U');
  jj = jj(jj > i);
  for j = jj
    x1 = P(1,i);
    y1 = P(2,i);
    x2 = P(1,j);
    y2 = P(2,j);
    x = linspace(x1+0.2*(x2-x1),x2+0.2*(x1-x2),numel(theta));
    x = x+mask*0.1*Radius(i).*cos(theta);
    y = linspace(y1+0.2*(y2-y1),y2+0.2*(y1-y2),numel(theta));
    y = y+mask*0.1*Radius(i).*sin(theta);
    plot(x,y,'k','LineWidth',1);
  end
end

drawnow;


%% a third time, defining domain and boundary voxels
figure(3), clf,
% these plot-commands makes the legend work...
adof = find(U);
idof = (N*(U ~= 0) > 0 & U == 0);
plot(P(1,adof),P(2,adof),'ks', ...
     P(1,idof),P(2,idof),'k.');

patch('Faces',R([adof; find(idof)],:), ...
      'Vertices',V,'FaceColor',[0.9 0.9 0.9]);
hold on,
axis([-1 1 -1 1]); axis square, axis off

ii = find(U == 1);
jj = find(U > 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'),'FaceAlpha',0.6);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'),'FaceAlpha',0.6);
axis([-0.25 1 -0.25 1]);

% ...here we plot again
adof = find(U);
idof = (N*(U ~= 0) > 0 & U == 0);
plot(P(1,adof),P(2,adof),'ks', ...
     P(1,idof),P(2,idof),'k.');
h = legend('\Omega_h','\partial \Omega_h');


%% a fourth time showing movement rates

% Set movement and pressure parameters
r_die = 1;      % rate of death -> pressure sink strength

% diffusive pressure rates
Drate1 = 1;    % into free matrix
Drate2 = 1;      % into already visited matrix
Drate3 = 0.1; % into already occupied voxel
Drate_ = [Drate1 Drate2; NaN Drate3];
% (hence Drate_(U(j)+1,VU(j)+1) = Drate_(2*VU(j)+U(j)+1) is the rate
% to move into voxel j)

% oversimplification given the unstructured mesh
% but OK for visualisation purposes 
gradquotient = 1; 

% boundary conditions: free matrix BC1, already visited matrix BC2 (0:
% no pressure, +: must overcome a pressure, -: cells get "sucked" out)
BC1 = 0;
BC2 = 0;

VU = (U > 0);   % visited voxels

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
    
% pressure Laplacian
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

% $$$ % pressure plot
% $$$ figure(3), clf,
% $$$ Pr_ = full(U); Pr_(Adof) = Pr;
% $$$ trisurf(T(1:3,:)',P(1,:),P(2,:),Pr_);

% visualize rates
figure(4), clf,
patch('Faces',R([adof; idof],:), ...      % boundary voxels
      'Vertices',V,'FaceColor',[0.9 0.9 0.9]);
hold on,
axis equal, axis off
ii = find(U == 1);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'), 'FaceAlpha', 0.6);
jj = find(U > 1);
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',graphics_color('vermillion'), 'FaceAlpha', 0.6);
axis([-0.25 1 -0.25 1]);

% visited boundary voxel or not
plot(P(1,idof1),P(2,idof1),'ko');
plot(P(1,idof2),P(2,idof2),'k.');

% compute all "moveb"-rates
x = zeros(1,0); y = zeros(1,0);
u = zeros(1,0); v = zeros(1,0);
position_scaling = 3; % scale factor to manually ensure arrows cross the edge
k_ = 0;
for i = bdof_m'
  k_ = k_+1;
  jj_ = find(N(i,Adof));
  keep = find(U(Adof(jj_)) == 0);
  jj_ = jj_(keep);
  for j_ = jj_
    if Pr(bdof_m_(k_)) > Pr(j_)
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient...
          *(Pr(bdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient...
          *(Pr(bdof_m_(k_))-Pr(j_))];
      x = [x 0.5*(P(1,i) + P(1,Adof(j_))) - position_scaling*u(end)];   
      y = [y 0.5*(P(2,i) + P(2,Adof(j_))) - position_scaling*v(end)];
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
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)...
          *gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)...
          *gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
      x = [x 0.5*(P(1,i) + P(1,Adof(j_))) - position_scaling*u(end)];
      y = [y 0.5*(P(2,i) + P(2,Adof(j_))) - position_scaling*v(end)];
    end
  end
end
quiver(x,y,u,v, 'k', 'LineWidth', 1.0,'AutoScaleFactor', 1.0);
drawnow;

%% Paper-specific formatting
figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 280 240]);
ax = gca;           % Tightens frame
ax.Position = ax.OuterPosition;
% uncomment to save:
%print -depsc 'figures/basic_modeling1.eps'

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 240]);
ax = gca;
ax.Position = ax.OuterPosition;
drawnow;
% uncomment to save:
%print -depsc 'figures/basic_modeling2.eps'

figure(3),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340 240]);
ax = gca;
ax.Position = ax.OuterPosition;
drawnow;
% uncomment to save:
%print -depsc 'figures/basic_modeling3.eps'

figure(4),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 280 240]);
ax = gca;
ax.Position = ax.OuterPosition;
drawnow;
% uncomment to save:
%print -depsc 'figures/basic_modeling4.eps'