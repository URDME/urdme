%STARTUP. Initializes URDME.
%
% B. Drawert   2010-06-07

function [] = urdme_startup() 

setenv('PATH',strcat( getenv('PATH'), ':/usr/local/bin' ));

%%% URDME initialization
if(isempty(getenv('URDME_ROOT')))
    fprintf(['URDME_ROOT environmental variable not set,' ... 
             'attempting to set it...\n']);
    [s,r] = system('urdme_init -r');
    uroot = strtok(r);
    if(isdir(uroot)),setenv('URDME_ROOT',uroot),end
    if(isempty(getenv('URDME_ROOT'))) %Only works from the examples dir
        [s,r] = system('../bin/urdme_init -r');
        uroot = strtok(r);
        if(isdir(uroot)),setenv('URDME_ROOT',uroot),end
    end
    
    if(isempty(getenv('URDME_ROOT'))) %Only works from the examples dir
        [s,r] = system('../urdme/bin/urdme_init -r');
        uroot = strtok(r);
        if(isdir(uroot)),setenv('URDME_ROOT',uroot),end
    end
    clear r s uroot;
    if(isempty(getenv('URDME_ROOT')))
        fprintf('URDME_ROOT environmental variable not set, exiting\n');
        return;
    else
        fprintf('URDME_ROOT set sucessfully.\n');
    end
end

% Add 'msrc' folder and all its subfolders to the path.
addpath(genpath(strcat(getenv('URDME_ROOT'),'/msrc')));
addpath(genpath(strcat(getenv('URDME_ROOT'),'/comsol')));

end