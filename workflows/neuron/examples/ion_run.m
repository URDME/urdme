%ION_RUN Rallpack 3 comparision.
%   Runs a sample simulation from Rallpack 3 with a time
%   discretization of 0.05 and a total of 100 voxels. This script also
%   gives an example on the use of the MODEL_BUILDER and the
%   NEURON_SOLVER functions, where tree structures are prepared for
%   discretization and solved using URDME.

% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-05-31

if ~exist('report','var')
  report = 1;
end

% first voxel reference voltage from Rallpack3
data = dlmread('examples/data/rallpack3_refsol1.txt');
V_pre = data(:,2)*1000;

% geometry
clear umod;
nVoxels = 100;
umod.private.model = model_builder(nVoxels);
% (default is 1mm long, 1um thick straight cable)

% Define input current (an nVoxel long array), I_inj = input current
% [nA] on the initial voxel after t [s].
%
% Replace (t >= 0) --> (mod(t,60) >= 45) to get a modulated current.
umod.private.I_inj = @(V,t,i)([0.1*(t >= 0); zeros(nVoxels-1,1)]);

% solve
dt = 0.05;
t_vec = 0:dt:250;
[V_out,I_mem,umod] = neuron_solver(t_vec,umod,'channel/',report);

% compute injected current
I_inj = zeros(size(t_vec));
for i = 1:numel(t_vec)-1
  Inj = umod.private.I_inj(V_out(i,1),t_vec(i),i);
  I_inj(i+1) = Inj(1);
end

% plots
if exist('plotting_off','var') && plotting_off, return; end
figure(1), clf,
plot(t_vec,V_out(:,1)); 
hold on
if exist('yyaxis','file')
  plot(t_vec,V_pre(1:end-1));
  ylabel('Membrane potential [mV]');
  yyaxis right
  plot(t_vec,I_inj,'LineStyle','--');  
  ylim([0 1.05*max(I_inj)]);
  ylabel('Injected current [nA]'); 
else
  [ax,h1,h2] = plotyy(t_vec,V_pre(1:end-1),t_vec,I_inj);
  set(h2,'LineStyle','--');
  set(get(ax(1),'ylabel'),'string','Membrane potential [mV]');
  set(get(ax(2),'ylabel'),'string','Injected current [nA]');
  set(ax(2),'ylim',[0 1.05*max(I_inj)]);
end
xlabel('Time [ms]');
legend('URDME','Rallpack 3','Initial current');

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[300 300 800 450]);
