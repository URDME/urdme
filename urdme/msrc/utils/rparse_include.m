function rparse_include(umod,includestr,where)
%RPARSE_INCLUDE Insert include-statement in file.
%   RPARSE_INCLUDE(UMOD,INCLUDESTR) inserts the string INCLUDESTR,
%   typically of the form '#include "foo.h"' into file UMOD.propensities
%   at a suitable place.
%
%   RPARSE_INCLUDE(FILENAME,INCLUDESTR) is also supported.

% Experimental syntax:
%   RPARSE_INCLUDE(UMOD,INCLUDESTR,WHERE) places INCLUDESTR according
%   to the string WHERE, with WHERE = 'include' (the default), or WHERE =
%   'propensity' understood as the propensity definition segment.

% S. Engblom 2022-10-19 (Revision, WHERE)
% S. Engblom 2019-12-03

if isstruct(umod)
  filename = umod.propensities;
else
  filename = umod;
end
if nargin < 3, where = 'include'; end

% read file
[fid,msg] = fopen(filename,'rt');
if fid == -1, error(msg); end  
F = char(fread(fid))';
if ~strncmp(F,'/* [Remove/modify this line not to overwrite this file] */',58)
  error('Will not overwrite existing file if not created by RPARSE.');
end
fclose(fid);

% insert includestr after...
switch where
 case 'include'
  s = '#include "report.h"';
 case 'propensity'
  s = '/* propensity definitions */';
 otherwise
  error('Placement specification not supported.');
end
ix = strfind(F,s);
if isempty(ix)
  error('Could not find place of insertion in file.');
end
F = [F(1:ix+length(s)) includestr F(ix+length(s):end)];

% write back
[fid,msg] = fopen(filename,'wt');
if fid == -1, error(msg); end
fprintf(fid,F);
fclose(fid);
