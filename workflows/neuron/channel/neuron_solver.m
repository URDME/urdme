function [V_out,I_mem,umod] = ...
    neuron_solver(tspan,umod,propwd,report)
%NEURON_SOLVER Split-step simulation of neuron.
%   [V_OUT,I_MEM,UMOD] = NEURON_SOLVER(TSPAN,UMOD,PROPWD,REPORT)
%   simulates a neuron in the time discretization TSPAN given model
%   data in UMOD. On return, V_OUT is the computed potential, I_MEM
%   the membrane current, and UMOD the updated model structure.
%
%   The model is described minimally by the fields UMOD.PRIVATE.model
%   and UMOD.PRIVATE.I_inj. The rest of the UMOD struct is then built
%   and returned.
%
%   Optional arguments include PROPWD, the working directory for
%   propensities, and REPORT, the report level of the solver.

% S. Engblom 2018-06-21 (Major revision)
% A. Senek 2017-05-31

% syntax
if nargin < 4
  report = 0;
  if nargin < 3
    propwd = '';
  end
end
if nargin < 2 || ~isfield(umod,'private') || ...
      ~isfield(umod.private,'model') || ~isfield(umod.private,'I_inj')
  error('Specify the neuron model in the field umod.private.');
end
model = umod.private.model;
nVoxels = model.nVoxels;

% default channels
if ~isfield(umod.private,'channels')
  umod.private.channels = squid_membrane(umod.private.model);
end
channels = umod.private.channels;

% C-N matrices
if ~isfield(umod.private,'model_matrx')
  c = channels.lchannel.G.*channels.lchannel.Erev;
  S1 = -(sum(channels.G_axial,2)+channels.G_m+channels.lchannel.G);
  S2 = sparse(channels.G_axial);
  C = sparse(1:nVoxels,1:nVoxels,channels.C_m);

  umod.private.model_matrx.c = c;
  umod.private.model_matrx.S1 = S1;
  umod.private.model_matrx.S2 = S2;
  umod.private.model_matrx.C = C;
else
  c = umod.private.model_matrx.c;
  S1 = umod.private.model_matrx.S1;
  S2 = umod.private.model_matrx.S2;
  C = umod.private.model_matrx.C;
end

% URDME set-up
try
  % validate umod: potentially compile and parse, but no solve
  umod = urdme(umod,'solve',0);
  if umod.compile
    unix(['mv mexssa.' mexext ' mexssa_neuron_solver.' mexext]);
    umod.solver = 'ssa_neuron_solver';
  end
  umod.solve = 1;
  umod.compile = 0;
  umod.parse = 0;
  % (note: under this syntax, umod.u0 is used and set_u0 bypassed)
catch
  % otherwise build umod from scratch
  umod = HHKNa_CHANNEL(umod);

  nSpecies = size(umod.N,1);
  umod.D = sparse(nSpecies*nVoxels,nSpecies*nVoxels);
  umod.vol = ones(1,nVoxels);
  umod.sd = ones(1,nVoxels);

  % calculate the number of ion-channels based on surface area
  ion_num = tprod(channels.density,model.surfacearea,[3 2],[1 3]);
  umod.u0 = set_u0(0,ion_num,propwd);
  umod.tspan = [0 1];
  umod.ldata = zeros(1,size(channels.Erest,1));

  % URDME parse and compile:
  umod = urdme(umod,'solver','ssa','solve',0);
  % avoid name clash with synaptic_solver:
  unix(['mv mexssa.' mexext ' mexssa_neuron_solver.' mexext]);
  umod.solver = 'ssa_neuron_solver';
  % solve next time:
  umod.solve = 1;
  umod.compile = 0;
  umod.parse = 0;
end

% open states
open_States = [5 13];
nOpenStates = numel(open_States);
num_open = zeros(nVoxels,nOpenStates);

% ODE stepping
Vm = umod.ldata'+channels.Erest;
V_out = zeros(numel(tspan),nVoxels);
V_out(1,:) = Vm';
I_mem = zeros(numel(tspan),nVoxels);

for i = 1:numel(tspan)-1
  % URDME simulation
  dt = tspan(i+1)-tspan(i);
  umod.tspan = tspan([i i+1]);
  umod.seed = randi(intmax('uint32'));
  umod = urdme(umod);

  % open states at end of interval
  xx = reshape(umod.U(:,2),[],nVoxels);
  umod.u0 = xx;
  num_open = xx(open_States,:)';

  % update matrices according to open channels
  c_ = c+num_open*(channels.G_single.*channels.Erev);
  S1_ = S1-num_open*channels.G_single;

  % C-N step
  t = tspan(i); % have not yet taken the step
  Iin_tmp = umod.private.I_inj(V_out(i,1),t,i);
  % (alternative: the average of I_inj() at tspan(i) and tspan(i+1)?)
  A = C/(dt/2)-(sparse(1:nVoxels,1:nVoxels,S1_)+S2);
  B = C/(dt/2)*Vm+(Iin_tmp+c_);
  Vm = 2*(A\B)-Vm;
  V_out(i+1,:) = Vm';
  I_mem(i+1,:) = Iin_tmp+c_+S1_.*Vm+S2*Vm;

  % prepare for next round
  umod.ldata = (Vm-channels.Erest)';

  % report progress
  if mod(i,round(numel(tspan)/10)) == 1 && report
    percent = floor(i/numel(tspan)*100);
    str = ['[' num2str(percent) '%%]-'];
    fprintf(str);
  end
end
if report, disp('[100%]'); end
