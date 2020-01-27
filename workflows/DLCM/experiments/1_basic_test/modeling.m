% Modeling explanation.
%
%   This script draws an illustrative schematic picture explaining the
%   numerical model.

% S. Engblom 2017-12-26 (revision)
% S. Engblom 2017-08-29 (minor revision)
% S. Engblom 2017-01-30

load experiments/petri_mesh/petri_geom % load [G,B]

% create mesh
[P,E,T] = initmesh(G,'hmax',0.35);
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
    for scale = [0.8 0.15];
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

% a second time, but with "connecting springs"
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

% a third time, defining domain and boundary voxels
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

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
%print -dpdf figures/basic_modeling1.pdf

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
drawnow;
% uncomment to save:
%print -dpdf figures/basic_modeling2.pdf

figure(3),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
drawnow;
% uncomment to save:
%print -dpdf figures/basic_modeling3.pdf
