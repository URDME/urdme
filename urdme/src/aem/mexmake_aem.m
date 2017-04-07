function mexmake_aem(propensity_file,~)
%MEXMAKE_AEM Makefile for MEXAEM.
%   MEXMAKE_AEM(P) Makes the AEM-solver with propensity source file P,
%   given as a relative path.

% S. Engblom 2017-02-17 (Major revision, URDME 1.3, Comsol 5)
% J. Cullhed 2008-06-18 (make)

if nargin > 1, error('AEM does not accept make arguments.'); end

% global defines, if any
define = [];

% path = location of this make
path = mfilename('fullpath');
path = path(1:end-11); % m-e-x-m-a-k-e-_-a-e-m is 11 chars

% determine path to propensity
if nargin > 0 && ~isempty(propensity_file)
  % propensity_file is a relative path
  propensity_source = [pwd '/' propensity_file];
else
  % can also compile mexnsm without propensity_file, using inline
  % propensities only
  propensity_source = [path '../propensities.c']';
end

% include and source directories
include = {['-I' path] ['-I' path '../../include']};
link =    {['-L' path] ['-L' path '../']};
source = {[path 'mexaem.c'] ...
          propensity_source ...
          [path 'aem.c'] ...
          [path '../binheap.c'] ...
          [path '../inline.c'] ...
          [path '../report.c']};
define = [define '-DMALLOC\(n\)=mxMalloc\(n\) -DFREE\(p\)=mxFree\(p\)'];

% mex extension
mx = mexext;

% platforms (edit here)
if strcmp(mx,'mexa64')
  cc = 'CC=gcc';
  cflags = ['CFLAGS=-fPIC -fno-omit-frame-pointer -std=c99 -O3 ' ...
            '-D_GNU_SOURCE -pthread -fexceptions '];
  
  mex('-silent','-largeArrayDims',cc,[cflags define], ...
      include{:},link{:},source{:});
elseif strcmp(mx,'mexmaci64')
  cflags = 'CFLAGS= -std=c99 ';
  mex('-silent','-largeArrayDims',[cflags define], ...
      include{:},link{:},source{:});
else
  error(['Platform not yet supported. Your MEX file extension is ' ...
         mx '. Please edit mexmake_aem.m to allow for this extension.']);
end
