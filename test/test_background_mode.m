% Integation test for background mode solver execution. 
%
% This test does not test the correctness of the NSM solver, 
% rather it tests that background mode execution works.  
%
% A. Hellander, 12/12/12


% Add mincde example to path

% Load previously stored input file
load('mincde_short_run.mat');

% Extend the mesh
umod = comsol2urdme(umod);

wd = pwd;
cd ../examples/mincde/
umod = mincde(umod);
% Make sure that it does not run for very long. 
umod.tspan = 0:20;
outfile = 'mincde_bg_run_out.mat';

% Passing an empty file-handle argument needs to work
umod = urdme(umod,[],'mode','bg','outfile',outfile);

while ~file_exists('mincde_bg_run_out.mat')
    sleep(1)
end

% Add solution back to the struct
umod = urdme_addsol(umod,outfile);

% Clean up
system(['rm ' outfile]);
cd wd;
