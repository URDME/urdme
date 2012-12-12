% URDME_SAVESOL Saves a trajectory stored in 'filename' from the FEM-struct. 
%
% B. Drawert, 2011-10-18. 
%

function fem = urdme_savesol(fem, filename)

%% Make sure that the urdme field is present
%if ~isfield(fem,'urdme')
%   error('No urdme field in fem-struct.'); 
%end

% Get trajectory from the model. 
U = fem.U;
tspan = fem.tspan;
save(filename,'U','tspan');


