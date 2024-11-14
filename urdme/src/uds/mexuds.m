function U = mexuds(mexhash,tspan,u0,D,N,G,vol,ldata,gdata, ...
                    data_time,ldata_time,gdata_time, ...
                    sd,reportl,seed,K,I,S,solverargs)

% (for help, type 'help uds')

% S. Engblom 2024-05-13 (data_time, ldata_time, gdata_time)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2019-11-17 (Revision, Nreplicas)
% S. Engblom 2017-02-24 (mexuds)
% S. Engblom 2016-01-07 (rt_solve, extension to linear transport models)
% S. Engblom 2007-04-05 (rrsolve)

% default solver options
optdef.odesolv = @ode23s;
optdef.odeopts = odeset('RelTol',1e-4,'AbsTol',1e-1);
optdef.mexname = 'mexuds';
optdef.report = 0;
optdef.jacobian = 0;
optdef.finish = '';

% input options
try
  opts = struct(solverargs{:});
catch
  error('Could not create Matlab struct from solver arguments.');
end
opts.report = reportl;

% merge defaults with actual inputs
fn = fieldnames(opts);
for i = 1:length(fn)
  if ~isfield(optdef,fn{i})
      error(sprintf('Unrecognized property ''%s''.',fn{i}));
    end
  optdef.(fn{i}) = opts.(fn{i});
end
opts = optdef;

if opts.report
  opts.odeopts.OutputFcn = @l_report;
end
if opts.jacobian
  opts.odeopts.Jacobian = @l_jacobian;
  NN = kron(eye(numel(vol)),N);
else
  % unused in this case:
  NN = N;
end

% create function handles
mexrhs = str2func([opts.mexname '_mexrhs']);
mexjac = str2func([opts.mexname '_mexjac']);

% solve
U = zeros(numel(u0(:,:,1)),numel(tspan),size(u0,3));
for k = 1:size(u0,3)
  l_report(0,[k size(u0,3)],'init_replica');
  [foo,U_] = opts.odesolv(@l_rhs,tspan,u0(:,:,k),opts.odeopts,mexhash, ...
                          N,NN,G,D,vol,ldata,gdata, ...
                          data_time,ldata_time,gdata_time, ...
                          sd,K,I,S, ...
                          mexrhs,mexjac);
  % fix for singular behaviour
  if numel(tspan) == 2
    U_ = U_([1 end],:);
  end
  U(:,:,k) = U_';
end

% play sound, if any
if ~isempty(opts.finish)
  load(opts.finish);
  sound(y,Fs);
end

%--------------------------------------------------------------------------
function dy = l_rhs(t,y,mexhash,N,NN,G,D,vol,ldata,gdata, ...
                    data_time,ldata_time,gdata_time, ...
                    sd,K,I,S, ...
                    mexrhs,mexjac)
%L_RHS Reaction-transport equations.
%  DY = L_RHS(...) returns the rate DY for the reaction-transport
%  model.

% if there is a table of times given, find the relevant place in
% (ldata_time,gdata_time)
if numel(data_time) > 1
  [~,ix] = histc(t,[data_time(:).' inf]);
  ldata_time = ldata_time(:,:,ix);
  gdata_time = gdata_time(:,ix);
end
dy = N*reshape(mexrhs(mexhash,t,y,size(N,2),vol,ldata,gdata, ...
                      ldata_time,gdata_time,sd,K,I,S), ...
               size(N,2),numel(vol));
dy = dy(:)+D*y;

%--------------------------------------------------------------------------
function J = l_jacobian(t,y,mexhash,N,NN,G,D,vol,ldata,gdata, ...
                        data_time,ldata_time,gdata_time, ...
                        sd,K,I,S, ...
                        mexrhs,mexjac)
%L_JACOBIAN Jacobian of L_RHS.
%  J = L_JACOBIAN(...) returns the Jacobian of L_RHS for the same
%  arguments.
if numel(data_time) > 1
  [~,ix] = histc(t,[data_time(:).' inf]);
  ldata_time = ldata_time(:,:,ix);
  gdata_time = gdata_time(:,ix);
end
J = NN*mexjac(mexhash,t,y,size(N,2),G,vol,ldata,gdata, ...
              ldata_time,gdata_time,sd,K,I,S)+D;

%--------------------------------------------------------------------------
function status = l_report(t,y,s,varargin)
%L_REPORT Simple reporter. Adapted from stenglib, see www.stenglib.org.

persistent T0 Tend percent t0 hwait kreplica Nreplicas;

if isempty(s) && ~isempty(hwait)
  % main use
  now = round(((kreplica-1)/Nreplicas+(t(end)-T0)/(Tend-T0)/Nreplicas)*100);
  if now > percent
    percent = now;
    if isempty(t0) || percent < eps
      waitbar(percent/100,hwait);
    end
  end
elseif strcmp(s,'init_replica')
  % init of each replica (must be performed before use)
  kreplica = y(1);
  Nreplicas = y(2);
elseif strcmp(s,'init')
  % main init of reporter (must be performed before use)
  T0 = t(1);
  Tend = t(end);
  percent = 0;
  t0 = [];
  try, close(hwait); catch, end
  hwait = waitbar(0.0,'Solution progress');
else
  % 'none' and 'done' clears the reporter
  Nreplicas = [];
  kreplica = [];
  try, close(hwait); catch, end
  hwait = [];
  t0 = [];
end
status = 0;

%--------------------------------------------------------------------------
