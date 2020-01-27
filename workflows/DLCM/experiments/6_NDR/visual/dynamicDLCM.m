%DYNAMICDLCM Visualizations of DYNAMIC_RUN.

% S. Engblom 2018-01-25

% uncomment in blocks to produce graphics

% $$$ load data/dynamic_growth
% $$$ 
% $$$ j = 0;
% $$$ for i = [100 500 2000 size(Etime,2)]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$ 
% $$$   U = Usave{1}+sum(Event(:,1:i),2);
% $$$   patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
% $$$   hold on,
% $$$   axis([-1 1 -1 1]); axis square, axis off
% $$$   ii = find(U == 1);
% $$$   patch('Faces',R(ii,:),'Vertices',V, ...
% $$$         'FaceColor',graphics_color('bluish green'));
% $$$   jj = find(U > 1);
% $$$   patch('Faces',R(jj,:),'Vertices',V, ...
% $$$         'FaceColor',graphics_color('vermillion'));
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamic_growth%d',j));
% $$$ end
% $$$ close all
% $$$ 
% $$$ return;

% $$$ load data/dynamicSSA1
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicSSA1_%d',j));
% $$$ end

% $$$ load data/dynamicSSA2
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicSSA2_%d',j));
% $$$ end

% $$$ load data/dynamicSSA6
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicSSA6_%d',j));
% $$$ end
% $$$ 
% $$$ return;

% $$$ load data/dynamicRDME1
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicRDME1_%d',j));
% $$$ end

% $$$ load data/dynamicRDME2
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicRDME2_%d',j));
% $$$ end

% $$$ load data/dynamicSSA6
% $$$ 
% $$$ mx = full(max(Dsave{end}(:)));
% $$$ j = 0;
% $$$ for i = [4 8 12 51]
% $$$   j = j+1;
% $$$   figure(j), clf,
% $$$   NDRplot(Dsave{i},Usave{i},V,R,mx);
% $$$ end
% $$$ 
% $$$ for j = 1:4
% $$$   figure(j)
% $$$   set(gcf,'PaperPositionMode','auto');
% $$$   set(gcf,'Position',[100 100 250 250]);
% $$$   axis equal, axis off
% $$$   % uncomment to save:
% $$$ % $$$   print('-depsc',sprintf('../../figures/dynamicRDME6_%d',j));
% $$$ end

return;
