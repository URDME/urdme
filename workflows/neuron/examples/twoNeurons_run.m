%TWONEURONS_RUN Complete coupling of two neurons.
%   Uses a reference solution from Rallpack 3 to propagate a potential
%   across a synapse into a neuron, and across a second synapse and
%   finally a second neuron.

% S. Engblom 2019-11-29 (Minor revision)
% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-05-31

if ~exist('report','var')
  report = 1;
end

% (1) synaptic solver using as input a Rallpack 3 potential
data = dlmread('examples/data/rallpack3_refsol2.txt');
t_vec = data(:,1).*1000; % [ms]
V_pre = data(:,2).*1000; % [mV]
if report, disp('Synaptic solver...'); end
[~,I_syn] = synaptic_solver([],t_vec,V_pre',report);
if report, disp('...done.'); end

% (2) neuron solver using I_syn as injected current
clear umod;
nVoxels = 30;
umod.private.model = model_builder(nVoxels,100,1);
umod.private.I_inj = @(V,t,i)((I_syn(i,:)+0.00) ...  % leakage current
                              .*[1; zeros(nVoxels-1,1)]);
if report, disp('Neuron solver...'); end
[V_out,~,umod] = neuron_solver(t_vec,umod,'channel/',report);
if report, disp('...done.'); end

% compute injected current
I_inj = zeros(size(t_vec));
for i = 1:numel(t_vec)-1
  Inj = umod.private.I_inj(V_out(i,1),t_vec(i),i);
  I_inj(i+1) = Inj(1);
end

if exist('plotting_off','var') && plotting_off, return; end

% frequency plot
T = 0.250; % [ms]
L = 5000;  % # datapoints
Fs = T/L;  % sampling freq

Y_pre = fft(V_pre);
P_temp = abs(Y_pre/length(Y_pre));
P_pre = P_temp(1:L/2+1);
P_pre(2:end-1) = 2*P_pre(2:end-1);

Y_out = fft(V_out(:,end));
P_temp = abs(Y_out/length(Y_out));
P_out = P_temp(1:L/2+1);
P_out(2:end-1) = 2*P_out(2:end-1);

f = Fs*(0:(L/2))/L;
figure(1), clf,
plot(f,P_pre)
title('Fourier spectra of input/output signals')
xlabel('Frequency');
ylabel('Amplitude');
hold on
plot(f,P_out);
legend('Input signal','Output signal');
xlim([0 2e-6])

% voltage plot
figure(2), clf,
subplot(2,1,1);
plot(t_vec,V_pre(:,end));
hold on
plot(t_vec,V_out(:,end));
legend('First neuron','Second neuron');
title('Membrane potential');
xlabel('Time [ms]');
ylabel('Voltage [mV]');

subplot(2,1,2);
plot(t_vec,I_inj);
title('Synaptic current')
xlabel('Time [ms]');
ylabel('Current [nA]');
