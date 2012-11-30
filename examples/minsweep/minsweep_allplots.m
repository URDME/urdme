%MINSWEEP_ALLPLOTS Superscript for all minsweep-related plots.
% S. Engblom 2011-11-21

load results/info % get xsep
x = cell(1,0);
y = cell(1,0);
p = cell(1,0);
for i = 1:numel(xsep)
  [xx,pp] = minsweep_plot(i);
  x = [x xx];
  y = [y repmat(xsep(i),1,numel(xx))];
  p = [p pp];
end

% illustrative of the case with no separation:
minsweep_plot(1);

% visualization of the point where reliable oscillations are no longer
% sustained:
figure(3),
for i = 1:numel(x)
  plot3(x{i},y{i},p{i},'k.--'); hold on
end
xlabel('x [\mu m]');
ylabel('sweep: separation [\mu m]');
zlabel('MinD_m [Relative concentration, temporal average]');
axis square
grid on