%BASIC_TEST_DRIVER Driver function to BASIC_TEST.

% S. Engblom 2017-08-29 (revision)
% S. Engblom 2016-12-25 (finalized physics - again!)
% S. Engblom 2016-12-09 (hexagonal mesh)
% S. Engblom 2016-11-14

clear all;
Nruns = 100;
report_progress = false;

% Cartesian
mesh_type = 1;
for irun = 1:Nruns
  irun
  rng(irun);
  basic_test;
  if irun == 1
    Ustats = Usave{end};
    P2Astats = P2A;
    P2Astats_var = P2A.^2;
  else
    Ustats = Ustats+Usave{end};
    P2Astats = P2Astats+P2A;
    P2Astats_var = P2Astats_var+P2A.^2;
  end
end
Ustats = full(Ustats/Nruns);
P2Astats = P2Astats/Nruns;
P2Astats_var = P2Astats_var/(Nruns-1)-Nruns/(Nruns-1)*P2Astats.^2;
V1 = V; R1 = R;

% hexagonal
mesh_type = 2;
for irun = 1:Nruns
  irun
  rng(irun);
  basic_test;
  if irun == 1
    Vstats = Usave{end};
    P2Astats_hex = P2A;
    P2Astats_hex_var = P2A.^2;
  else
    Vstats = Vstats+Usave{end};
    P2Astats_hex = P2Astats_hex+P2A;
    P2Astats_hex_var = P2Astats_hex_var+P2A.^2;
  end
end
Vstats = full(Vstats/Nruns);
P2Astats_hex = P2Astats_hex/Nruns;
P2Astats_hex_var = P2Astats_hex_var/(Nruns-1)-Nruns/(Nruns-1)*P2Astats_hex.^2;
V2 = V; R2 = R;

% uncomment to save:
save('experiments/1_basic_test/basic_test','Ustats','Vstats', ...
     'V1','R1','V2','R2','tspan','P2Astats','P2Astats_var', ...
     'P2Astats_hex','P2Astats_hex_var');

return;

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
colorbar;
% uncomment to save:
% $$$ print -depsc figures/basic_test.eps

figure(2),
clf,
h = patch('Faces',R2,'Vertices',V2,'FaceVertexCData',Vstats, ...
          'FaceColor','flat','EdgeColor','none');
hold on,
axis([-sqrt(3)/2 sqrt(3)/2 -1 1]); axis equal
i = find(~any(isnan(R2),2),1,'first');
a = polyarea(V2(R2(i,:),1),V2(R2(i,:),2));
r = sqrt(sum(Vstats(:))*a/pi);
z = r*exp(2i*pi*linspace(0,1,40));
plot(z,'k');
set(gca,'xtick',[],'ytick',[],'Visible','off');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
% $$$ print -depsc figures/basic_test_hex.eps

figure(3), clf,
errorbar(tspan,P2Astats,sqrt(P2Astats_var),'b.-');
hold on,
errorbar(tspan,P2Astats_hex,sqrt(P2Astats_hex_var),'r.-');
xlabel('time');
ylabel('P2A');
legend('Cartesian','Hexagonal');
axis tight
ylim([1 1.3]);
set(gca,'ytick',1:0.1:1.3,'xtick',[0 100]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 260]);
% uncomment to save:
% $$$ print -depsc figures/basic_test_P2A.eps