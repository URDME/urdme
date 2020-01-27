%NEURON_RUN Example run of realistic model neuron.
%   Runs a simulation using a pre-processed tree.

% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-05-31

if ~exist('report','var')
  report = 1;
end

% model and geometry parameters (#compartments in cable, [um])
nVoxels = 200;
% contains "tree"; voxel #2072 is middle:
load('examples/data/model.mat');
% find all dendrite voxels:
dends = find(tree.R == 3);
% find dendrite starter voxels:
starterV = find(sum(tree.dA(dends,dends),1) == 0);
% connect over the root to propagate directly from dendrites:
tree.dA(1183,2) = 1;
clear umod;
umod.private.model = tree;

% define input current
I_vec = zeros(tree.nVoxels,1);
I_vec(starterV) = 1;
umod.private.I_inj = @(V,t,i)(0.1*(t >= 0)*I_vec);

% solve
dt = 0.05;
t_vec = 0:dt:50;
[V_out,I_mem,umod] = neuron_solver(t_vec,umod,'channel/',report);

% compute injected current
I_inj = zeros(size(t_vec));
for i = 1:numel(t_vec)-1
  Inj = umod.private.I_inj(V_out(i,1),t_vec(i),i);
  I_inj(i+1) = Inj(1);
end

% plots
if exist('plotting_off','var') && plotting_off, return; end
if exist('yyaxis','file')
  yyaxis right
  plot(t_vec,V_out(:,1)); 
  ylabel('Membrane potential [mV]');
  hold on
  yyaxis left
  plot(t_vec,I_inj,'LineStyle','--');  
  ylabel('Injected current [nA]'); 
else
  [ax,h1,h2] = plotyy(t_vec,V_out(:,1),t_vec,I_inj);
  set(get(ax(1),'ylabel'),'string','Membrane potential [mV]');
  set(get(ax(2),'ylabel'),'string','Injected current [nA]');
  set(h2,'LineStyle','--');
end
xlabel('Time [ms]');

return;

% produce animation (requires TREES Toolbox)
make_gif(tree,t_vec,V_out,filename);
