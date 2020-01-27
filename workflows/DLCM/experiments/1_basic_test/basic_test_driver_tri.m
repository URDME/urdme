%BASIC_TEST_DRIVER_TRI Driver function to BASIC_TEST_TRI.

% S. Engblom 2017-01-14 (_tri)
% S. Engblom 2016-12-25 (finalized physics - again!)
% S. Engblom 2016-12-09 (hexagonal mesh)
% S. Engblom 2016-11-14

clear all;
Nruns = 100;
report_progress = false;

% unstructured triangulation
geom = 3;
for irun = 1:Nruns
  irun
  rng(irun);
  basic_test_tri;
  if irun == 1
    Ustats = Usave{end};
  else
    Ustats = Ustats+Usave{end}; 
  end
  reuse = true;
end
Ustats = full(Ustats/Nruns);
V1 = V; R1 = R;
clear reuse;

% uncomment to save:
% $$$ save('experiments/1_basic_test/basic_test_tri','Ustats', ...
% $$$      'V1','R1');

% visualization
figure(1), 
clf,
h = patch('Faces',R1,'Vertices',V1,'FaceVertexCData',Ustats, ...
          'FaceColor','flat','EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis equal
i = find(~any(isnan(R1),2),1,'first');
a = polyarea(V1(R1(i,:),1),V1(R1(i,:),2));
r = sqrt(sum(Ustats(:))*a/pi);
z = r*exp(2i*pi*linspace(0,1,40));
plot(z,'k');
set(gca,'xtick',[],'ytick',[],'Visible','off');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
% $$$ print -depsc figures/basic_test_tri.eps
