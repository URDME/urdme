function mexmake_uds(propensity_file,varargin)
%MEXMAKE_UDS Makefile for MEXRHS.
%   MEXMAKE_UDS(P) Makes the UDS-solver with propensity source file P,
%   given as a relative path.
%
%   MEXMAKE_UDS(P,...) accepts additional arguments as
%   property/value-pairs.
%
%   Property     Value/{Default}     Description
%   -----------------------------------------------------------------
%   define       string              Defines on the mex-format 
%                                    '-D<define>', see MEX.
%   source       string              Source filename(s).
%   include      string              Include path.
%   link         string              Linker path.

% S. Engblom 2020-02-21 (Revision, mexjac)
% S. Engblom 2019-12-03 (Revision, make arguments)
% S. Engblom 2019-11-06 (Revision, now using URDMEstate_t)
% S. Engblom 2017-02-24

% global defines, if any
define = [];

% path = location of this make
path = mfilename('fullpath');
path = path(1:end-11); % m-e-x-m-a-k-e-_-u-d-s is 11 chars

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

% include and source directories
include = {['-I' path] ['-I' path '../../include']};
if ~isempty(opts.include)
  include = [include ['-I' opts.include]];
end
link =    {['-L' path] ['-L' path '../']};
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
          '-DMALLOC\(n\)=mxMalloc\(n\) -DFREE\(p\)=mxFree\(p\)' ...
          opts.define];

% mex extension
mx = mexext;

% platforms (edit here)
if strcmp(mx,'mexa64')
  cc = 'CC=gcc';
  cflags = ['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
            '-D_GNU_SOURCE -pthread -fexceptions '];
  
  source{1} = sourceMEX{1};
  mex('-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
  source{1} = sourceMEX{2};
  mex('-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
elseif strcmp(mx,'mexmaci64')
  cflags = 'CFLAGS= -std=c99 -mmacosx-version-min=10.15 ';
  source{1} = sourceMEX{1};
  mex('-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
  source{1} = sourceMEX{2};
  mex('-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
else
  error(['Platform not yet supported. Your MEX file extension is ' ...
         mx '. Please edit mexmake_uds.m to allow for this extension.']);
end
