function [umod,I_syn] = synaptic_solver(umod,tspan,V_pre,report)
%SYNAPTIC_SOLVER Simulation of current across synaptic cleft.
%   [UMOD,I_SYN] = SYNAPTIC_SOLVER(UMOD,TSPAN,V_PRE,REPORT) calculates
%   the synaptic current in the time discretization TSPAN given the
%   potential V_PRE in the dendrite. If V_PRE is a matrix of size
%   M-by-N, then V_PRE(i,j) is the potential a time TSPAN(j) and
%   synaptic connection i.
%
%   The model is described by UMOD which can be empty on input in
%   which case the UMOD struct is built and returned.
%
%   REPORT is an optional argument, the report level of the
%   solver. The default is REPORT = 1.

% S. Engblom 2019-11-29 (Revision, new syntax)
% S. Engblom 2018-06-21 (Major revision)
% A. Senek 2017-05-31

% syntax
if nargin < 4
  report = 1;
end

nVoxels = size(V_pre,1); 
open_indx = [4 11];
num_tot = 1000;

g_max_AMPA = 1;     % [nS]
g_max_GABAA = -0.1; % [nS]
g_max = [g_max_AMPA g_max_GABAA];

% calcium release
L_max = 1.54; % [mM]
Vp = 2;       % [mV]
Kp = 5;       % [mV]
umod.private.conc_func = @(V)L_max./(1+exp(-(V-Vp)./Kp));

try
  % validate umod: potentially compile and parse, but no solve
  umod = urdme(umod,'solve',0,'modelname','synaptic_solver');
  if umod.compile
    %unix(['mv mexssa.' mexext ' mexssa_synaptic_solver.' mexext]);
    umod.solver = 'ssa_synaptic_solver';
  end
  umod.solve = 1;
  umod.compile = 0;
  umod.parse = 0;
  % (note: under this syntax there is no URDME "burn in")
catch
  % build umod
  umod = AMPA_NMDA_SYNAPSE(umod);

  nSpecies = size(umod.N,1);
  umod.D = sparse(nVoxels*nSpecies,nVoxels*nSpecies);
  umod.vol = ones(1,nVoxels);
  umod.sd = ones(1,nVoxels);

  % initial data for propensity function
  umod.ldata = zeros(2,nVoxels);
  umod.ldata(1,:) = V_pre(1);                         % "v"
  umod.ldata(2,:) = umod.private.conc_func(V_pre(1)); % "T"

  % set the initial distribution of synaptic molecules (needs a
  % reference for the number num_tot of ions available in synapse)
  umod.u0 = zeros(nSpecies,nVoxels);
  umod.u0(1,:) = num_tot;
  umod.u0(7,:) = 0.7*num_tot*ones(1,nVoxels);
  umod.u0(12,:) = 0.3*num_tot*ones(1,nVoxels);
    
  % URDME synapse simulation: "burn in"
  umod.tspan = [0 10];
  umod.seed = randi(intmax('uint32'));
  umod = urdme(umod,'solver','ssa','modelname','synaptic_solver');
  umod.u0 = reshape(umod.U(:,2),nSpecies,nVoxels);

  % avoid name clash with neuron_solver:
  %unix(['mv mexssa.' mexext ' mexssa_synaptic_solver.' mexext]);
  umod.solver = 'ssa_synaptic_solver';
  umod.compile = 0;
  umod.parse = 0;
end

% consecutive solve
I_syn = zeros(numel(tspan),nVoxels);
for i = 1:numel(tspan)-1
  umod.ldata(1,:) = V_pre(:,i);                         % "v"
  umod.ldata(2,:) = umod.private.conc_func(V_pre(:,i)); % "T"
  umod.tspan = tspan([i i+1]);
  umod.seed = randi(intmax('uint32'));
  umod = urdme(umod);

  xx = reshape(umod.U(:,2),[],nVoxels);  
  I_tmp = (g_max*(xx(open_indx,:)./num_tot))';
  I_syn(i+1,:) = I_tmp;

  % prepare for next round
  umod.u0 = xx;

  % report progress
  if mod(i,round(numel(tspan)/10)) == 1 && report
    percent = floor(i/numel(tspan)*100);
    str = ['[' num2str(percent) '%%]-'];
    fprintf(str);
  end
end
if(report), disp('[100%]'); end
