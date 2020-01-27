%CIRCLENEURONS_RUN Complete coupling of neurons.
%   Uses a reference solution from Rallpack 3 to propagate a potential
%   across 9 straight neurons, which are connected by synapses. This
%   result can then be used to create a spatially propagating model
%   using Comsol.

% S. Engblom 2019-11-29 (Minor revision)
% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-06-13

if ~exist('report','var')
  report = 1;
end

load('data/circle_tree.mat');
prop_relpath_chnl = 'channel/';

% break the connection between segments
input_nodes = [16 29 44 57 71 84 97 111 2];
pre_syn_nodes = [14 27 42 55 69 82 95 109 123];
tree.dA(input_nodes,:) = 0;

% setup neuron channel model
clear umod_chnl;
umod_chnl.private.model = tree;
nVoxels = tree.nVoxels;
input_vec = zeros(1,nVoxels);
umod_chnl.private.I_inj = @(V,t,i) ((t > 0).*input_vec)';

% create all necessary matrices
t_setup = [0 0.01];
[~,~,umod_chnl] = neuron_solver(t_setup,umod_chnl,prop_relpath_chnl);

% setup synaptic channels
bogus_Vm = zeros(length(input_nodes),1)-75;
[umod_syn,~] = synaptic_solver([],t_setup,bogus_Vm,0); 

% ODE stepping
dt = 0.05;
t_vec = 0:dt:25;
V_out = zeros(numel(t_vec)-1,nVoxels);
V_out(1,:) = -75;
Iin = zeros(numel(t_vec)-1,length(input_nodes));
I_out = zeros(numel(t_vec)-1,nVoxels);

for i = 1:numel(t_vec)-1
  t_tmp = t_vec([i i+1]);
  % channel simulation
  [V_temp,I_mem,umod_chnl] = neuron_solver(t_tmp,umod_chnl,prop_relpath_chnl);

  % synaptic simulation, set the pre-synaptic voltage (2 voxels before
  % the post-synaptic voxel)
  V_pre = V_temp(2,pre_syn_nodes)';
  [umod_syn,I_syn] = synaptic_solver(umod_syn,t_tmp,V_pre,0);

  % update synaptic current
  inj_vec = input_vec;
  inj_vec(input_nodes) = I_syn(2,:);
  umod_chnl.private.I_inj = @(V,t,i)((t > 0).*inj_vec)'+ ...
      (0.5*(t < 0.2).*[1 zeros(1,nVoxels-1)])';

  Iin_tmp = umod_chnl.private.I_inj(V_out(i,1),t_tmp(2),i);
  Iin(i+1,:) = Iin_tmp(input_nodes);
  V_out(i+1,:) = V_temp(2,:);
  I_out(i+1,:) = I_mem(2,:);
  I_out(i+1,input_nodes-1) = I_syn(2,:);

  % report progress
  if mod(i,round(numel(t_vec)/10)) == 1 && report
    percent = floor(i/numel(t_vec)*100);
    str = ['[' num2str(percent) '%%]-'];
    fprintf(str);
  end
end
if(report), disp('[100%]'); end

%  voltage plot
if exist('plotting_off','var') && plotting_off, return; end
figure(1), clf
subplot(2,1,1);
plot(t_vec,V_out(:,pre_syn_nodes));
title('Neuron voltage in the neurons (9 in total)');
xlabel('Time [ms]');
ylabel('Voltage [mV]');

subplot(2,1,2);
plot(t_vec,Iin);
title('Synaptic current from respective neuron');
xlabel('Time [ms]');
ylabel('Current [nA]');

return;

% save the solution for COMSOL use
I = [t_vec; I_out']';
save('spatial/I.mat','I');

% animation
load('data/tree.mat');
make_gif(umod.private.model,t_vec,V_out,'circle_fire.gif');
