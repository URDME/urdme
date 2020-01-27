function U = mexuds(tspan,u0,D,N,G,vol,ldata,gdata,sd,reportl,seed,K,I,S,solverargs)

% (for help, type 'help uds')

% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2019-11-17 (Revision, Nreplicas)
% S. Engblom 2017-02-24 (mexuds)
% S. Engblom 2016-01-07 (rt_solve, extension to linear transport models)
% S. Engblom 2007-04-05 (rrsolve)

% default solver options
optdef.odesolv = @ode23s;
optdef.odeopts = odeset('RelTol',1e-4,'AbsTol',1e-1);
optdef.report = 0;
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
  optdef = setfield(optdef,fn{i},getfield(opts,fn{i}));
end
opts = optdef;
if opts.report
  opts.odeopts.OutputFcn = @l_report;
end

% solve
U = zeros(numel(u0(:,:,1)),numel(tspan),size(u0,3));
for k = 1:size(u0,3)
  l_report(0,[k size(u0,3)],'init_replica');
  [foo,U_] = opts.odesolv(@l_rhs,tspan,u0(:,:,k),opts.odeopts,N,D,vol, ...
                          ldata,gdata,sd,K,I,S);
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
function dy = l_rhs(t,y,N,D,vol,ldata,gdata,sd,K,I,S)
%L_RHS Reaction-transport equations.
%  DY = L_RHS(T,Y,N,D,VOL,...) returns the rate DY for the
%  reaction-transport model at time T. The system is described by the
%  reactions N and the transport rates D in voxel volumes VOL.

dy = reshape(N*mexrhs(t,y,size(N,2),vol,ldata,gdata,sd,K,I,S),[],1)+D*y;

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
