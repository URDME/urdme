% RDME2MAT. Write field in fem.urdme to file. 
%
% A. Hellander, 2010-04-27.
%

function rdme2mat(umod,filename)

urdme_validate(umod); 
%don't save the fem-model
if ~isfield(umod,'test')
  umod=rmfield(umod,'comsol');
end
save(filename, '-struct', 'umod','-v6');