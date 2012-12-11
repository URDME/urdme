function [umod,outputfile] = urdme(fem,varargin)
%URDME Interface to spatial stochastic simulation algorithms.
%   UMOD = URDME(COMSOL,FILE) generates a single trajectory using the model
%   .m-file FILE as the specification of the model. It is assumed the
%   propensity functions are found in the model .c-file with the same
%   name.
%
%   For example, UMOD = URDME(COMSOL,'mymodel') uses the model .m-file
%   'mymodel.m' by executing UMOD = mymodel(COMSOL) and also compiles the
%   model propensity .c-file 'mymodel.c'.
%
%   UMOD = URDME(COMSOL,@FUN,'propensities','mymodel') uses the function
%   handle FUN instead and the propensity file 'mymodel.c' for the
%   reactions.
%
%   UMOD = URDME(COMSOL,FILE,...) passes additional options in the form of
%   property/value pairs, see the table below.
%
%   For example, after installing the DFSP-solver, you may use the
%   syntax
%     UMOD = URDME(COMSOL,FILE,'solver','dfsp')
%   to launch the solver named 'dfsp' in interactive mode. 
%
%   Solvers can generally be run in background mode by specifying the
%   property 'mode' to be 'bg', that is,
%     [UMOD,OUTFILE] = URDME(COMSOL,FILE,'mode','bg')
%   will launch the default NSM solver in background mode, returning
%   control to the Matlab-prompt immediately. See the example below.
%
%   Property     Value/{Default}        Description
%   -----------------------------------------------------------------------
%   solver       {'nsm'} | ...          Name of solver.
%   mode         {'interactive'} | 'bg' Mode of operation.
%   verbose      {0}, 1, 2              Level of Matlab layer report.
%   report        0, {1}, 2             Level of solver report.   
%   seed         uint32                 Random number seed.
%   propensities character array        Name of propensity .c-file.
%   outfile      character array        Name of output file
%   solvopts     character array        Solver options
%   delete_inputfile    0,{1}                  Delete the input file. 
%   delete_outputfile   0,{1}                  Delete the output file. 
%
%   The default NSM solver supports three report levels, 0 (no
%   report), 1 (progress report), 2 (progress report + event count).
%
%   The URDME matlab interface supports three levels of verbosity. 
%   0 (silent), 1 (intermediary) and 2 (comprehensive). Level 1 can be 
%   of interests to an advanced user, while the highest level
%   is only recommened as a tool to debug the interface. 
%
%   Example:
%     [fem,out] = urdme(fem,file,'mode','bg'); % background mode
%     % wait until done
%     load(out);
%
%     % Comsol 3.5-syntax:
%     fem = rdme2fem(fem,U);
%
%     % Comsol 4.x-syntax:
%     fem = rdme2mod(fem,model,U);
%     % (where model is the Comsol Java-object)
%
%   See also MOD2RDME, FEM2RDME, RDME2MOD, RDME2FEM, RDME2MAT.

% B. Drawert, A. Hellander 2012-10-2 (Revision, inline propensities)
% P. Bauer and S. Engblom 2012-04-04 (Revision, cleanup)
% V. Gerdin 2012-02-01 (Revision, Comsol 4.2)
% B. Drawert and A. Hellander 2010-06-07 (Revision, background mode)
% J. Cullhed 2008-06-18

% make sure everything is initialized
urdme_startup;

% valid options (supported by NSM-solver) 
optdef = struct('report',1,'seed',[],'mode','interactive', ...
                'solver','nsm','outfile',[],'propensities',[], ...
                'solvopts',[],'verbose',1,'delete_inputfile',1, ...
                'delete_outputfile',1);

% parse property/value pairs
if nargin > 2 
  try
    % URDME 1.1 syntax
    if iscell(varargin{2})
      arr=varargin{2};
      opts = struct(arr{:});  
    % URDME 1.2 syntax
    else
      opts = struct(varargin{2:end});
    end
  catch
    error('Options must be passed as property/value pairs.');
  end
    
  % merge fields
  fn = fieldnames(opts);
  for i = 1:length(fn)
      % In general, we canot restrict what options are passed in to the
      % solvers. New solvers may need new options.  Also all DFSP options must be supported.
      % This is a consequence of the overall design of the interface layer.
      % To warn the user, we instead issue a warning if verbose > 0. 
      if ~isfield(optdef,lower(fn{i}))
            warning(['Property name ''' fn{i} ''' is not a valid option for the core nsm solver.']);
      end
      optdef = setfield(optdef,lower(fn{i}),getfield(opts,fn{i}));
  end
end

opts = optdef;

% check for urdme struct
if ~isfield(fem,'test')
  if nargin < 2
    error(['No .urdme field in fem-struct and no model ' ...
           'function handle specified.']);
  else
    umod=comsol2urdme(fem,opts.verbose);
  end
else
  % externally defined umod structure (for tests)
  umod=fem;
end

% We don't impose any restrictions on report levels (all integers) since
% add-on solvers may suport a highel report level than 2. You should
% write report level logic tests using ">", not "==".  
if isnumeric(opts.report) && opts.report==int32(opts.report)
    umod.report = opts.report;
else
    fprintf('report=%g\n',opts.report);
    error('Unsupported report level (must be integer).');
end

umod.seed = opts.seed;

if strcmp(opts.mode,'bg')
    md = '&';
elseif ~strcmp(opts.mode,'interactive')
    error('Unsupported mode.');
else
    md = '';
end

% execute .m-model file, 3.5/4.x-style (empty propensity allowed for
% the purpose of convenient testing)
if nargin > 1 && ~isempty(varargin{1})
  if ~isa(varargin{1},'function_handle')
    if ~exist([varargin{1} '.m'],'file')
      error(['Could not find model file ' ...
             varargin{1} ...
             '.m in current directory or path.']);
    end
  end
  % execute model file
  umod = feval(varargin{1},umod);
end

% name of model file
if ischar(varargin{1})
  model_name = varargin{1};
else
  model_name = '';
end
% If no model file is given, assume that it is named the
% same as the current working directory. 
if(isempty(model_name))
    [t,r]=strtok(pwd,'/');
    while ~isempty(r)
        [t,r]=strtok(r,'/');
    end
    model_name = t;
end


% check that fem.urdme contains all data structures required by the
% NSM solver and that they make sense
urdme_validate(umod);

% precompiled propensities
if isempty(opts.propensities) 
    if(~isempty(model_name) && file_exists([model_name,'.c']))
        opts.propensities = strcat(model_name,'.c');
    end
elseif(~file_exists(opts.propensities) && file_exists(strcat(opts.propensities,'.c')))
        opts.propensities = strcat(opts.propensities,'.c');
end

% add 'msrc' directory to path, if using alternate solver path
if strcmp(getenv('URDME_SOLVER_PATH'),'')
    rest = getenv('URDME_SOLVER_PATH');
    while(~isempty(rest))
        [tok,rest]=strtok(rest,':');
        if(isdir([tok,'/src/',opts.solver]))
            orig_path = path();
            path(orig_path,[tok,'/msrc/']);
            if(opts.verbose>=2)
                fprintf('adding %s to path\n',[tok,'/msrc/']);
            end
            break;
        end
    end
end

% initialize solvers using an optional initialization script.
s = strcat('urdme_init_',char(opts.solver));
if exist(s)
    if opts.verbose>=2
      fprintf('executing %s\n',s);
    end
    umod = eval(strcat(s,'(umod,varargin)'));
% elseif ~strcmp(getenv('URDME_SOLVER_PATH'),'')
%     rest = getenv('URDME_SOLVER_PATH');
%     while(~isempty(rest))
%         [tok,rest]=strtok(rest,':');
%         if(isdir([tok,'/src/',opts.solver]))
%             fid=fopen([tok,'/msrc/urdme_init_',opts.solver,'.m']);
%             if(fid~=-1)
%                 orig_path = path();
%                 path(orig_path,[tok,'/msrc/']);
%                 if(exist(s))
%                     if(opts.verbose>=2)
%                         fprintf('\texecuting %s in %s\n',s,[tok,'/msrc/']);
%                     end
%                     umod = eval(strcat(s,'(umod,varargin)'));
%                 else
%                     error(sprintf('file %s is found, but exits(%s)==0\n',s,s));
%                 end
%                 path(orig_path);
%             end
%             break;
%         end
%     end
end

% get safe temporary files for input/output
[foo,inputfile] = system('mktemp -t urdmemodel.XXXXXXXXXX');
inputfile = strcat(strtrim(inputfile), '.mat');
if opts.verbose>1
  fprintf('inputfile=%s\n',inputfile);
end

% if output file not supplied by the user, get a temporary file
if isempty(opts.outfile)
  [foo,opts.outfile] = system('mktemp -t urdmeresults.XXXXXXXXXX');
  opts.outfile = strcat(strtrim(opts.outfile), '.mat');
end
if opts.verbose>1
  fprintf('outputfile=%s\n',opts.outfile);
end

% compile the solver
if opts.verbose > 0
  fprintf('Compiling solver...\n');
end

umod = urdme_compile(umod,model_name,opts.propensities,opts.solver,opts.verbose);

%fem = urdme_compile(fem,opts.propensities,opts.solver,opts.solvopts);
if opts.verbose ~= 0
  fprintf('   ...done.\n\n');
end

% Serialize model to temporary input-file.
if(opts.verbose>1)
    fprintf('Writing temporary input file.\n',inputfile);
end
rdme2mat(umod,inputfile);
% solve!
if opts.verbose ~= 0
  fprintf('Starting simulation...\n');
end
cmd = ['./.urdme/' model_name '.' opts.solver ' ' inputfile ' ' ...
       opts.outfile ' ' md];
   
if(opts.verbose>1)
    fprintf('cmd=%s\n',cmd)
end  
system(cmd);

% if not in background mode, load the result file and add the solution
% to the field fem.urdme
if ~strcmp(opts.mode,'bg')
    
  if ~exist(opts.outfile,'file')
    system(['rm ' inputfile]);
    error('Solver did not finish correctly.');
  else 
    data=load(opts.outfile);
    
    % if a test, a manually created mesh or model with sd=1 has been executed,
    % only add U to umod, otherwise add solution to comsol file as well. 
    %if isfield(umod,'test') || length(unique(umod.sd))==1
      umod.U = data.U;
      if isfield(data,'tspan') 
          umod.tspan = data.tspan;
      else
          umod.tspan = 0:size(data.U,2)-1;
      end
      umod = urdme2comsol(umod,data.U,umod.tspan,opts.verbose);
  end


  % clean temporary input/output files.
  if(opts.delete_inputfile)
      if(opts.verbose>2)
          fprintf('rm %s\n',inputfile);
      end
      system(['rm ' inputfile]);
  end
  if(opts.delete_outputfile)
      if(opts.verbose>2)
          fprintf('rm %s\n',opts.outfile);
      end
      system(['rm ' opts.outfile]);
  end
  
end

outputfile = opts.outfile;
