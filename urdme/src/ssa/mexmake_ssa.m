function mexmake_ssa(propensity_file,varargin)
%MEXMAKE_SSA Makefile for MEXSSA.
%   MEXMAKE_SSA(P) Makes the SSA-solver with propensity source file P,
%   given as a relative path.
%
%   MEXMAKE_SSA(P,...) accepts additional arguments as
%   property/value-pairs.
%
%   Property     Value/{Default}     Description
%   -----------------------------------------------------------------
%   OpenMP       Boolean {false}     Turns OpenMP-compilation
%                                    on/off.
%   define       string              Defines on the mex-format 
%                                    '-D<define>', see MEX.
%   source       string              Source filename(s).
%   include      string              Include path.
%   link         string              Linker path.

% S. Engblom 2019-11-29 (Revision, make arguments)
% S. Engblom 2017-02-22

% global defines, if any
define = [];

% path = location of this make
path = mfilename('fullpath');
path = path(1:end-11); % m-e-x-m-a-k-e-_-s-s-a is 11 chars

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
optdef.openmp = false;
optdef.define = '';
optdef.include = '';
optdef.link = '';
optdef.source = '';

% merge defaults with actual inputs
if nargin > 1
  opts = struct(varargin{:});
  fn = fieldnames(opts);
  for i = 1:length(fn)
    optdef = setfield(optdef,fn{i},getfield(opts,fn{i}));
  end
end
opts = optdef;

% OpenMP
if opts.openmp
  omp_link = '-lgomp';
  omp_cflags = '-fopenmp ';
else
  omp_link = '';
  omp_cflags = '';
end

% include and source directories
include = {['-I' path] ['-I' path '../../include']};
if ~isempty(opts.include)
  include = [include ['-I' opts.include]];
end
link = {omp_link ['-L' path] ['-L' path '../']};
if ~isempty(opts.link)
  link = [link ['-L' opts.link]];
end
source = {[path 'mexssa.c'] ...
          propensity_source ...
          [path 'ssa.c'] ...
          [path '../inline.c'] ...
          [path '../report.c']};
if ~isempty(opts.source)
  if ischar(opts.source)
    source = [source {opts.source}];
  else
    source = [source opts.source];
  end
end
define = [define '-DMALLOC\(n\)=mxMalloc\(n\) -DFREE\(p\)=mxFree\(p\) ' ...
          opts.define];

% mex extension
mx = mexext;

% platforms (edit here)
if strcmp(mx,'mexa64')
  cc = 'CC=gcc';
  cflags = ['CFLAGS= -fPIC ' omp_cflags ...
            '-fno-omit-frame-pointer -std=c99 -O3 ' ...
            '-D_GNU_SOURCE -pthread -fexceptions '];
  
  mex('-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
elseif strcmp(mx,'mexmaci64') || strcmp(mx,'mexmaca64') 
  if opts.openmp
    warning('OpenMP not supported on this platform.');
    link = link(2:end); % current fix: simply remove omp_link
  end
  cflags = 'CFLAGS= -std=c99 -mmacosx-version-min=10.15 ';
  mex('-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
else
  error(['Platform not yet supported. Your MEX file extension is ' ...
         mx '. Please edit mexmake_ssa.m to allow for this extension.']);
end
