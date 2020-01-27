function rparse_include(filename,includestr)
%RPARSE_INCLUDE Insert #include-statement in file.
%   RPARSE_INCLUDE(FILENAME,INCLUDESTR) inserts the string INCLUDESTR,
%   typically of the form '#include "foo.h"' into file FILENAME on a
%   suitable place.
  
% S. Engblom 2019-12-03
    
% read file
[fid,msg] = fopen(filename,'rt');
if fid == -1, error(msg); end  
F = char(fread(fid))';
if ~strncmp(F,'/* [Remove/modify this line not to overwrite this file] */',58)
  error('Will not overwrite existing file if not created by RPARSE.');
end
fclose(fid);

% construct modified version
s = '#include "report.h"'; % insert includestr after this
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
