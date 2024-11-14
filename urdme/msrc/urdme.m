function umod = urdme(umod,varargin)
%URDME Interface to spatial stochastic simulation algorithms.
%   UMOD = URDME(UMOD,...) generates a single trajectory for the
%   reaction-diffusion/transport model defined in the URDME struct
%   UMOD. The arguments to URDME are passed via property/value pairs,
%   see the table below.
%
%   For example,
%     UMOD = URDME(UMOD,'solver','mysolver','propensities','mymodel')
%   parses the URDME struct UMOD, compiles the solver 'mysolver',
%   assumed to have the make-file 'make_mysolver.m', together with the
%   propensity source file 'mymodel.c' for the reactions. The
%   resulting executable is then called using the arguments provided
%   in UMOD.
%
%   Comsol Java objects are conventionally stored in UMOD.comsol and
%   PDE Toolbox structures in UMOD.pde. After simulating, the
%   resulting object may be loaded back into Comsol (or PDE Toolbox)
%   for visualization and postprocessing. See URDME2COMSOL and
%   URDME2PDE.
%
%   Property        Value/{Default}         Description
%   -----------------------------------------------------------------------
%   solver          string {'nsm'}          Name of solver
%   propensities    string                  Path to propensity .c-file
%   report          {0}, 1, 2, 3            Level of report
%   solve           {1} | 0                 Solve on/off
%   compile         {1} | 0                 Compile on/off
%   parse           {1} | 0                 Parse on/off
%   seed            double {0}              Solver random number seed
%   modelname       string {''}             Name of model
%
%   The URDME Matlab interface supports three levels of report: 0
%   (silent), 1 (intermediate), 2 (comprehensive), 3 (as 2, but also
%   allowing for an early exit).
%
%   Turn compilation off when you are solving the same model several
%   times and after the first call to URDME. You may similarly
%   explicitly turn solving/parsing on/off. For example:
%     % first call, compile and solve mymodel.c:
%     UMOD = URDME(UMOD,'propensities','mymodel','seed',1234);
%     % second call, new seed to parse but no compilation:
%     UMOD = URDME(UMOD,'seed',2345,'compile',0);
%     % a sequence of calls:
%     UMOD.parse = 0; % don't parse the URDME struct
%     for seed = 1:10
%       UMOD.seed = seed; % set seed directly
%       UMOD = URDME(UMOD);
%     end
%
%   The property modelname gives control of the name of the executable
%   file, and is sent as an argument to the makefile, e.g.,
%     mex -outfile mexname ...
%   where mexname is constructed as 
%     mexname = mex<solver>[_<modelname>][_<propensities>]
%   Consider to set modelname if the same URDME struct is parsed and
%   compiled multiple times so as to provide for unique executable
%   names.
%
%   The URDME solver sequence takes the generic form
%     % if UMOD.compile:
%     UMOD.mexexec = make_<UMOD.solver>(UMOD.propensities,UMOD.makeargs);
%   Note that mexname and mexhash are considered standard input make
%   arguments, i.e., are always present.
%     % if UMOD.solve:
%     UMOD.U =  FEVAL(UMOD.mexexec,UMOD.mexhash, ...
%                     UMOD.tspan,UMOD.u0, ...
%                     UMOD.D,UMOD.N,UMOD.G, ...
%                     UMOD.vol,UMOD.ldata,UMOD.gdata,UMOD.sd, ...
%                     UMOD.report,UMOD.seed,UMOD.solverargs);
%
%   See also COMSOL2URDME, URDME2COMSOL, URDME_VALIDATE, NSM.

%   The URDME-struct has the following fields:
%   ------------------------------------------
%
%   Required fields before call:
%
%   tspan                Time vector
%   u0                   Initial state
%   D                    Diffusion matrix
%   N                    Stoichiometric matrix
%   G                    Dependency graph
%   vol                  Voxel volumes
%   sd                   Subdomain numbers
%
%   Required, usually passed as options to URDME:
%
%   solver               Solver
%   propensities         Propensity source file
%   report               Solver feedback
%   solve                Solve on/off
%   compile              Compilation on/off
%   parse                Parsing on/off
%   seed                 Random seed
%
%   Optional, empty understood when left out:
%
%   modelname            Name of model
%   inline_propensities  Inline propensity structure
%   ldata                Local data vector
%   gdata                Global data vector
%   data_time            Nodes for time-dependent data
%   ldata_time           Local time-dependent data vector
%   gdata_time           Global time-dependent data vector
%   solverargs           Arguments to solver, property/value cell-vector
%   makeargs             Arguments passed to makefile
%
%   Optional:
%
%   U                    Latest stored solution
%   comsol               Comsol Java object
%   pde                  PDE Toolbox structure
%   private              Field to store anything extra
%
%   Handled by URDME:
%
%   mexname              Name of executable mex-file
%                        {mex<solver>[_<modelname>][_<propensities>]}
%   mexhash              Hashkey of model (input to both makefile
%                        and solver)
%                        {automatic}
%   mexexec              Function handle of executable mex-file
%                        {automatic}

% S. Engblom 2024-05-04 (Revision, modelname/mexname/mexhash/mexexec)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% P. Bauer, S. Engblom 2012-04-04 (Revision, cleanup)
% V. Gerdin 2012-02-01 (Revision, Comsol 4.2)
% B. Drawert, A. Hellander 2010-06-07 (Revision, background mode)
% J. Cullhed 2008-06-18

%% (0) Initial parse of property/value pairs.

% the logic of this if-statement is: (1) either options have been
% given, or (2) the umod-struct is incomplete, or (3) the umod-struct
% itself instructs us to parse it
if nargin > 1 || ~isfield(umod,'parse') || umod.parse
  % default options
  optdef = struct('solver','nsm', ...
                  'propensities','', ...
                  'inline_propensities', ...
                  struct('K',[],'I',[],'S',[]), ...
                  'report',0, ...
                  'solve',1, ...
                  'compile',1, ...
                  'parse',1, ...
                  'seed',0, ...
                  'modelname','', ...
                  'ldata',[], ...
                  'gdata',[], ...
                  'data_time',[], ... 
                  'ldata_time',[], ...
                  'gdata_time',[], ...
                  'solverargs',{{}}, ...
                  'makeargs',{{}}, ...
                  'mexname','', ...
                  'mexhash',[], ...
                  'mexexec',[], ...
                  'U',[], ...
                  'comsol',[], ...
                  'pde',[], ...
                  'private',[]);
  % required fields, no meaningful defaults:
  req = struct('tspan',[],'u0',[],'D',[],'N',[],'G',[],'vol',[],'sd',[]);

  % input options
  try
    opts = struct(varargin{:});
  catch
    error('Could not create Matlab struct from input property/value pairs.');
  end

  % (1) merge options: input opts --> umod, opts takes precedence
  fn = fieldnames(opts);
  for i = 1:length(fn)
    umod.(fn{i}) = opts.(fn{i});
  end

  % (2) merge options: umod --> optdef, umod takes precedence
  fn = fieldnames(umod);
  for i = 1:length(fn)
    if ~isfield(optdef,fn{i}) && ~isfield(req,fn{i})
      error(sprintf('Unrecognized property ''%s''.',fn{i}));
    end
    if ~strcmp(fn{i},'inline_propensities')
      optdef.(fn{i}) = umod.(fn{i});
    else
      % the field inline_propensities is a struct itself and it is
      % convenient to have its empty default fields propagate into any
      % partially built struct
      fn_ = fieldnames(umod.inline_propensities);
      for j = 1:length(fn_)
       if ~isfield(optdef.inline_propensities,fn_{j})
         error(sprintf(['Unrecognized property ' ...
                        '''inline_propensities.%s''.'],fn_{j}));
       end
        optdef.inline_propensities.(fn_{j}) = umod.inline_propensities.(fn_{j});
      end
    end
  end
  % done!
  umod = optdef;
end

%% (1) Parse the input: possibly add/reformat missing default fields
% and then validate the resulting URDME struct
if umod.parse
  if isempty(umod.ldata)
    umod.ldata = zeros(0,numel(umod.vol));
  end
  if isempty(umod.data_time)
    umod.data_time = umod.tspan(1);
  end
  if isempty(umod.ldata_time)
    umod.ldata_time = zeros(0,numel(umod.vol),numel(umod.data_time));
  end
  if isempty(umod.gdata_time)
    umod.gdata_time = zeros(0,numel(umod.data_time));
  end
  % propensities, if any
  if ~isempty(umod.propensities) && ~any(umod.propensities == '.')
    umod.propensities = [umod.propensities '.c'];
  end
  % determine mexname from inputs {solver, modelname, propensities}:
  % always prepend solver name...
  if isempty(umod.modelname)
    umod.mexname = ['mex' umod.solver];
  else
    umod.mexname = ['mex' umod.solver '_' umod.modelname];    
  end
  % ...and always append name of propensity file (if any)
  if ~isempty(umod.propensities)
    name = umod.propensities;
    if any(name == '/')
      name = name(find(name == '/',1,'last')+1:end);
    end
    name = name(1:find(name == 'c',1,'last')-2);
    umod.mexname = [umod.mexname '_' name];
  end
  % in summary: mexname = mex<solver>[_<modelname>][_<propensities>]
  if isempty(umod.mexhash)
    % will do something more involved here in future releases:
    umod.mexhash = 0x6327d1c4; % (CRC-32 of 'URDME')
  end
  urdme_validate(umod);
else
  l_info(umod.report,2,'Parsing turned off.\n');
end

%% (2) Compile the solver.
if umod.compile
  mexmake = ['mexmake_' umod.solver];
  if exist(mexmake) ~= 2
    warning('URDME:noMake', ...
            sprintf('It seems there is no make for solver = ''%s.''.', ...
                    umod.solver));
  end
  % mexname & mexhash are always passed in order to support
  % identification of executables:
  umod.mexexec = feval(mexmake,umod.propensities, ...
                       'mexname',umod.mexname,'mexhash',umod.mexhash, ...
                       umod.makeargs{:});
else
  l_info(umod.report,2,'Compilation turned off.\n');
end

%% (3) Solve!
if umod.solve
  if umod.report >= 2, solver_timer = tic; end 
  l_info(umod.report,1,'Starting simulation...\n');
  exist_mexexec = exist(func2str(umod.mexexec));
  solverargs = umod.solverargs;
  if exist_mexexec == 2
    % singular behavior for .m-file solvers (e.g., UDS) since the actual
    % solver is not compiled; we might need to be able to find other
    % compiled files (e.g., mexrhs/mexjac), so we add the mexname as a
    % solver argument in order to be able to do this:
    solverargs = [solverargs 'mexname' umod.mexname];
  end
  if exist_mexexec ~= 2 && exist_mexexec ~= 3
    warning('URDME:noMex', ...
            sprintf(['It seems there is no .mex/.m-executable ' ...
                     'for solver = ''%s''.\n'],umod.solver));
  end
  umod.U = feval(umod.mexexec, ...
                 umod.mexhash, ...
                 umod.tspan,umod.u0,umod.D,umod.N,umod.G, ...
                 umod.vol,umod.ldata,umod.gdata, ...
                 umod.data_time,umod.ldata_time,umod.gdata_time, ...
                 umod.sd, ...
                 umod.report,umod.seed, ...
                 umod.inline_propensities.K, ...
                 umod.inline_propensities.I, ...
                 umod.inline_propensities.S, ...
                 solverargs);
  l_info(umod.report,1,'   ...done.\n');
  if umod.report >= 2
    fprintf('Solver execution time = %gs.\n',toc(solver_timer));
  end
else
  l_info(umod.report,2,'Solving turned off.\n');
end

%-------------------------------------------------------------------------
function l_info(level,lim,msg)
%L_INFO Display information.
%   L_INFO(LEVEL,LIM,MSG) Displays the message MSG whenever LEVEL >=
%   LIM.

if level >= lim, fprintf(msg); end

%-------------------------------------------------------------------------
