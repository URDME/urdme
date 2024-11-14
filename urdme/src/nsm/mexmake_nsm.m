function mexexec = mexmake_nsm(propensity_file,varargin)
%MEXMAKE_NSM Makefile for MEXNSM.
%   MEXEXEC = MEXMAKE_NSM(P) Makes the NSM-solver with propensity
%   source file P, given as a relative path, and returns the
%   executable file as a function handle.
%
%   MEXEXEC = MEXMAKE_NSM(P,...) accepts additional arguments as
%   property/value-pairs.
%
%   Property     Value/{Default}     Description
%   -----------------------------------------------------------------
%   mexname      string {'mexnsm'}   Name of output mex-function, cf.
%                                    "mex -output <mexname>"
%   mexhash      scalar uint32 {0}   Model hashkey, cf.
%                                    "-DMEXHASH=<mexhash>"
%
%   define       string              Defines on the mex-format 
%                                    '-D<define>', see MEX.
%   source       string              Source filename(s).
%   include      string              Include path.
%   link         string              Linker path.
%
%   See also NSM, MEX.

% S. Engblom 2024-05-03 (Revision, mexname/mexhash/mexexec)
% S. Engblom 2017-02-16 (Major revision, URDME 1.3, Comsol 5)
% J. Cullhed 2008-06-18 (make)

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
optdef.mexname = 'mexnsm';
optdef.mexhash = 0x0;
optdef.define = '';
optdef.source = '';
optdef.include = '';
optdef.link = '';

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
link = {['-L' path] ['-L' path '../']};
source = {[path 'mexnsm.c'] ...
          propensity_source ...
          [path 'nsm.c'] ...
          [path '../binheap.c'] ...
          [path '../inline.c'] ...
          [path '../report.c']};
define = [define '-DMALLOC\(n\)=mxMalloc\(n\) -DFREE\(p\)=mxFree\(p\) ' ...
          sprintf('-DMEXHASH=%u ',opts.mexhash) ...
          opts.define];

% mex extension
mx = mexext;

% platforms (edit here)
if strcmp(mx,'mexa64')
  cc = 'CC=gcc';
  cflags = ['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
            '-D_GNU_SOURCE -pthread -fexceptions '];
  mex('-output',opts.mexname,'-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
elseif strcmp(mx,'mexmaci64') || strcmp(mx,'mexmaca64') 
  cflags = 'CFLAGS= \$CFLAGS -std=c99 -mmacosx-version-min=10.15 ';
  mex('-output',opts.mexname,'-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
else
  error(['Platform not yet supported. Your MEX file extension is ' ...
         mx '. Please edit mexmake_nsm.m to allow for this extension.']);
end
mexexec = str2func(opts.mexname);
