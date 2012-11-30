% Adds a trajectory stored in 'filename' to a corresponding FEM struct. 
%
% Deals with several special cases that primarily arises when one
% want to add a solution stored on file to an exisiting fem-struct that
% may have been read from file or is newly exported from Comsol. In both 
% those cases the fem-struct may be missing filds that leads rdme2fem to 
% fail. 
%
% In normal circumstances, fem = urdme_addsol(fem,filename)
% is doing the same thing as 
%
% load(filename)
% fem.urdme.U = U;
% fem.urdme.tspan = tspan;
% fem = rdme2fem(fem);
%
%
% A. Hellander, 2010-06-12. 
%
% 

function umod = urdme_addsol(umod, filename, verbose)

data=load(filename);

% Check for endended mesh. 
if ~isfield(umod,'xmesh')||~umod.xmesh.initialized
   % If the urdme field is present, meshextend will choke. 
   
   if isfield(umod.comsol,'mesh')&&isfield(umod.comsol,'equ');   
     
       % PB: can you change this to the new structure?
       % umod=fem.urdme
       % umod.comsol=fem (without .urdme field)
       
       urdmefield = fem.urdme;
       fem = rmfield(fem,'urdme');
       fem.xmesh = meshextend(fem);
       fem.urdme=urdmefield;
   else
       if isfield(fem,'mesh')&&isfield(fem,'equ')
          fem.xmesh=meshextend(fem);
       else % If not using Comsol. 
          umod.U = data.U;
          return;
       end
   end
end
s
% Make sure that the urdme field is present. It may not be if we are adding
% a solution to a fem-struct that we have loaded from file, or a 
% newly exported struct from Comsol. 

% PB: In the new structure, there is no fem.urdme field, so what is it
% really you want to check for?
if ~isfield(fem,'urdme')
   fem = fem2rdme(fem);
end

% Add trajectory to model. In addition to adding the concentation data
% to the comsol data structure, we also save the raw copy numbers in the
% struct as they are frequently needed by postproicessing routines. 
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
% rdme2fem does not work (and is not needed). Here we simply add the 
% U matrix to the struct. 
if numel(fem.urdme.sd) > 1
    umod = urdme2comsol(umod,data.U,verbose);
end