%LAYERS Visualize the computational layers.

% S. Engblom 2018-01-24

% RDME layer
[P,E,T,G] = RDMElayer;

figure(1), clf, h = pdegplot(G); hold on
set(h,'linewidth',2,'color','k')
triplot(T(1:3,:)',P(1,:),P(2,:),'k');
axis tight, axis equal, axis off
 
% dual mesh
[V,R] = mesh2dual(P,E,T,'voronoi');
figure(2), clf,
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9], ... 
      'EdgeColor','black');
axis tight, axis equal, axis off

% DLCMlayer
figure(3), clf,
patch('Faces',R,'Vertices',V, ...
      'FaceColor',graphics_color('bluish green'), ... 
      'EdgeColor',[0 0 0]);
hold on
h = pdegplot(G);
set(h,'linewidth',2,'color','k')
axis tight, axis equal, axis off
V_ = V;
R_ = R;
P_ = P;

figure(4), clf,
[P,E,T,~,V,R] = DLCMlayer(20);
ii = find(sqrt(P(1,:).^2+P(2,:).^2) < 0.4);
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8]);
patch('Faces',R(ii,:), ...
      'Vertices',V,'FaceColor',graphics_color('bluish green'), ...
      'EdgeColor',[0 0 0]);
axis([-0.5 0.5 -0.5 0.5]), axis equal, axis off

figure(5), clf,
% same again, but zoomed in
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8]);
patch('Faces',R(ii,:), ...
      'Vertices',V,'FaceColor',graphics_color('bluish green'), ...
      'EdgeColor',[0 0 0],'LineStyle',':');
axis equal, axis([-0.5 0 -0.5 0]), axis off
for jj = [129:131 148:150 168:171 187:190]
  patch('Faces',R_,'Vertices',tsum(V_*0.05,P(:,jj),[1 2],[2 3]), ...
        'FaceColor',graphics_color('bluish green'), ... 
        'EdgeColor','black');
end

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
% $$$ print -depsc ../../figures/RDMEtri.eps

figure(2),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
% $$$ print -depsc ../../figures/RDMEdual.eps

figure(3),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 345 240]);
% uncomment to save:
print -depsc ../../figures/RDMElayer.eps

figure(4),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 345 240]);
% uncomment to save:
print -depsc ../../figures/DLCMlayer.eps

figure(5),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 345 240]);
% uncomment to save:
print -depsc ../../figures/DLCM_RDME.eps