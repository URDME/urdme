%STATICSSA Visualization of STATIC_RUN.

% S. Engblom 2018-01-25

load data/staticSSA

j = 0;
mx = max(Dsave{end}(:));
for i = [2 6 11 51]
  j = j+1;
  figure(j), clf,
  NDRplot(Dsave{i},Usave{i},V,R,mx);
end

for j = 1:4
  figure(j)
  set(gcf,'PaperPositionMode','auto');
  set(gcf,'Position',[100 100 250 250]);
  axis equal, axis off
  % uncomment to save:
% $$$   print('-depsc',sprintf('../../figures/staticSSA%d',j));
end
