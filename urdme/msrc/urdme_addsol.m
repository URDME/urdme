% Add a trajectory stored in 'filename' to a corresponding UMOD struct. 
%
% Deals with special cases that primarily arises when one
% want to add a solution stored on file to an exisiting umod-struct that
% may have been read from file or is newly exported from Comsol. In both 
% those cases the umod.comsol field may be missing, which leads 
% comsol2urmde to fail. 
%
% In normal circumstances, fem = urdme_addsol(fem,filename)
% is doing the same thing as 
%
% load(filename)
% fem.U = U;
% fem.tspan = tspan;
% fem = urdme2comsol(fem);
%
%
% A. Hellander, 2010-06-12. 
%
% 

function umod = urdme_addsol(umod, filename, verbose)

if nargin < 3
    verbose = 0;
end

data=load(filename);

% Add trajectory to model. In addition to adding the concentation data
% to the Comsol data structure, we also save the raw copy numbers in the
% struct as they are frequently needed by postprocessing routines. 
if isfield(data,'U') 
    umod.U = data.U;
else
    error(sprintf('No trajectory found in %s\n',filename));
end

% If we are adding a solution file to a previously saved fem-struct, 
% for example if we have been running jobs in background more, 
% we also add the tspan-vector. If no such vector is present in the 
% output file, tspan is set to unit intervals. 
if isfield(data,'tspan') 
    umod.tspan = data.tspan;
elseif ~isfield(umod,'tspan')
    umod.tspan = 0:size(data.U,2)-1;
end

% In the special case of only one subdomain (a well mixed simulation), 
% urdme2comsol does not work (and is not needed). 
if numel(umod.sd)>1&&isfield(umod,'comsol')
    umod = urdme2comsol(umod,umod.U,umod.tspan,verbose);
end
