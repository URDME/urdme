%STARTUP Startup for URDME.

% E. Blom 2025-12-12 (dlcm startup)
% S. Engblom 2017-02-14

% link = location of this startup.m
link = mfilename('fullpath');
link = link(1:end-numel(mfilename)); % remove 'startup'
if ~strcmp(pwd,link(1:end-1))
  warning('URDME startup should run from the URDME-folder.');
end

% path to urdme/ folders
addpath(genpath([link 'urdme/']));
addpath(genpath([link 'workflows/']));
addpath(genpath([link 'examples/']));

% startup DLCM
cd workflows/DLCM/
startup
cd(link)
