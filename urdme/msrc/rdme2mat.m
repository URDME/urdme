% RDME2MAT. Write field in fem.urdme to file. 
%
% A. Hellander, 2010-04-27.
%

function rdme2mat(umod,filename)

urdme_validate(umod); 
% We do not save the comsol datastructure
if ~isfield(umod,'test') && isfield(umod,'comsol');
  umod=rmfield(umod,'comsol');
end
save(filename, '-struct', 'umod','-v6');