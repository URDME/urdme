%PROTRUSIONS Sample visualization of protrusion model.

% S. Engblom 2018-01-24 (minor revision)
% S. Engblom 2017-08-30 (minor revision)
% S. Engblom 2017-04-03

rng(1704); % repeatable results

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 10;

% build the hex-mesh
[P,E,T] = basic_mesh(2,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');
Nvoxels = size(P,2);

% radius taken to be max distance from all vertices to voxel center
Radius = R;
for j = 1:size(R,1)
  % all NaN's mapped to first vertex
  Radius(j,isnan(Radius(j,:))) = R(j,1);
end
Radius = sqrt(max(sum(tsum(reshape(V(Radius,:),Nvoxels,size(R,2),2),-P, ...
                           [1 3 2],[2 1]).^2,2),[],3));

figure(1), clf,
ii = ~any(isnan(R),2);
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.3 0.3 0.3]);
hold on,

ii = 34;
patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',[207 24 20]/255,'EdgeColor',[0.3 0.3 0.3]);

jj = 36;
patch('Faces',R(jj,:),'Vertices',V, ...
      'FaceColor',[207 24 20]/255,'EdgeColor',[0.3 0.3 0.3]);

contact(P(:,ii),Radius(ii),4*Radius(ii),pi/4, ...
        P(:,jj),Radius(jj),4*Radius(jj),-pi/4,pi/12,1);

kk = 57;
patch('Faces',R(kk,:),'Vertices',V, ...
      'FaceColor',[207 24 20]/255,'EdgeColor',[0.3 0.3 0.3]);
contact(P(:,kk),Radius(kk),4*Radius(kk),pi*0.3, ...
        P(:,jj),Radius(jj),4*Radius(jj),-pi/4,pi/12,2);

axis equal, axis tight, axis off
h1 = text(-0.53,-0.28,'A');
h2 = text(0.36,-0.27,'B');
h3 = text(0.25,0.28,'C');
h4 = text(0.64,0.28,'\theta');
h5 = text(0.64,0.65,'d\theta');

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 400 285]);
% uncomment to save:
% $$$ print -depsc ../../figures/protrusions.eps

return;

% to determine cell #
p = ginput(1);
[~,ii] = min(sum(tsum(P,-p,[1 2],[3 1]).^2))
