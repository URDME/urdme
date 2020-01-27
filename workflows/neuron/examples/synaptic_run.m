%SYNAPTIC_RUN Current propagation across synaptic cleft.
%   Uses a reference voltage from Rallpack 3 to propagate a current
%   across a synapse by AMPAa and NMDAa channels.

% S. Engblom 2019-11-29 (Minor revision)
% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-05-31

if ~exist('report','var')
  report = 1;
end

% end-voxel voltage
data = dlmread('data/destexhe_refsol.txt');
t_vec = data(201:1300,1).*1000;
V_pre = data(201:1300,2).*1000;  

V_pre = [V_pre(end).*ones(100,1); V_pre; V_pre(end).*ones(300,1)]';
t_vec = linspace(0,35,1500);
[umod,num_open] = synaptic_solver([],t_vec,V_pre,report);

% calculate injected current
g_max_AMPA = 0.4; % [nS]
g_max_GABAA = 0;  % [nS]
g_max = [g_max_AMPA g_max_GABAA];
I_syn = tsum(tprod(g_max,num_open,[4 3],[1 2]),2:3); % leakage current

% plots
if exist('plotting_off','var') && plotting_off, return; end
figure(1), clf
subplot(3,1,1);
plot(t_vec,V_pre);
title(['Input voltage ']);
xlabel('Time [ms]');
ylabel('Voltage [mV]');
ylim([-100 50])

subplot(3,1,2);
plot(t_vec,umod.private.conc_func(V_pre));
title('Transmitter concentration');
xlabel('Time [ms]');
ylabel('Concentration [mM]');

subplot(3,1,3);
plot(t_vec,-I_syn');
title('Synaptic current');
xlabel('Time [ms]');
ylabel('Current [nA]');
ylim([-0.15 0]);
