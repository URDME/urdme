function mexexec = mexmake_uds(propensity_file,varargin)
%MEXMAKE_UDS Makefile for MEXRHS.
%   MEXEXEC = MEXMAKE_UDS(P) Makes the UDS-solver with propensity
%   source file P, given as a relative path, and returns the
%   executable file as a function handle.
%
%   MEXEXEC = MEXMAKE_UDS(P,...) accepts additional arguments as
%   property/value-pairs.
%
%   Property     Value/{Default}     Description
%   -----------------------------------------------------------------
%   mexname      string {'mexuds'}   Name of output mex-function, cf.
%                                    "mex -output <mexname>"
%                (note: '_mexrhs/_mexjac' appended)
%   mexhash      scalar uint32 {0}   Model hashkey, cf.
%                                    "-DMEXHASH=<mexhash>"
%
%   define       string              Defines on the mex-format 
%                                    '-D<define>', see MEX.
%   source       string              Source filename(s).
%   include      string              Include path.
%   link         string              Linker path.
%
%   Use 'define' = '-DJAC_' to compile and line the interface to a
%   prepared Jacobian for compiled propensities (this is not required
%   for inline propensities).
%
%   See also UDS, MEX.

% S. Engblom 2024-05-03 (Revision, mexname/mexhash/mexexec)
% S. Engblom 2020-02-21 (Revision, mexjac)
% S. Engblom 2019-12-03 (Revision, make arguments)
% S. Engblom 2019-11-06 (Revision, now using URDMEstate_t)
% S. Engblom 2017-02-24

% global defines, if any
define = [];

% path = location of this make
path = mfilename('fullpath');
path = path(1:end-numel(mfilename));

% determine path to propensity
if nargin > 0 && ~isempty(propensity_file)
  % propensity_file is a relative path
  propensity_source = [pwd '/' propensity_file];
else
  % can also compile without propensity_file, using inline propensities
  % only
  propensity_source = [path '../propensities.c'];
end

% default options
optdef.mexname = 'mexuds';
optdef.mexhash = 0x0;
optdef.define = '';
optdef.include = '';
optdef.link = '';
optdef.source = '';

% merge defaults with actual inputs
if nargin > 1
  opts = struct(varargin{:});
  fn = fieldnames(opts);
  for i = 1:length(fn)
    try
      % trigger warning on missing field?
      f = optdef.(fn{i});
      % if not: accept the field
      optdef.(fn{i}) = opts.(fn{i});
    catch
      warning(sprintf('Make argument ''%s'' seems unsupported.',fn{i}));
    end
  end
end
opts = optdef;

% include and source directories
include = {['-I' path] ['-I' path '../../include']};
if ~isempty(opts.include)
  include = [include ['-I' opts.include]];
end
link = {['-L' path] ['-L' path '../']};
if ~isempty(opts.link)
  link = [link ['-L' opts.link]];
end
% compiles two mexFunction():
sourceMEX = {[path 'mexrhs.c'] ...
             [path 'mexjac.c']};
source = {[] ... % (place to insert one of the above)
          propensity_source ...
          [path '../inline.c']};
if ~isempty(opts.source)
  if ischar(opts.source)
    source = [source {opts.source}];
  else
    source = [source opts.source];
  end
end
define = [define '-DUDS_ ' ...
          '-DMALLOC\(n\)=mxMalloc\(n\) -DFREE\(p\)=mxFree\(p\) ' ...
          sprintf('-DMEXHASH=%u ',opts.mexhash) ...
          opts.define];

% mex extension
mx = mexext;

% platforms (edit here)
if strcmp(mx,'mexa64')
  cc = 'CC=gcc';
  cflags = ['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
            '-D_GNU_SOURCE -pthread -fexceptions '];
  source{1} = sourceMEX{1};
  mex('-output',[opts.mexname '_mexrhs'], ...
      '-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
  source{1} = sourceMEX{2};
  mex('-output',[opts.mexname '_mexjac'], ...
      '-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
elseif strcmp(mx,'mexmaci64') || strcmp(mx,'mexmaca64') 
  cflags = 'CFLAGS= \$CFLAGS -std=c99 -mmacosx-version-min=10.15 ';
  source{1} = sourceMEX{1};
  mex('-output',[opts.mexname '_mexrhs'], ...
      '-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
  source{1} = sourceMEX{2};
  mex('-output',[opts.mexname '_mexjac'], ...
      '-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
else
  error(['Platform not yet supported. Your MEX file extension is ' ...
         mx '. Please edit mexmake_uds.m to allow for this extension.']);
end
mexexec = str2func('mexuds');