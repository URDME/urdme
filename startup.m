%STARTUP Startup for URDME.

% S. Engblom 2017-02-14

% link = location of this startup.m
link = mfilename('fullpath');
link = link(1:end-numel(mfilename)); % remove 'startup'

% path to urdme/ folders
addpath(genpath([link 'urdme/']));
addpath(genpath([link 'workflows/']));
addpath(genpath([link 'examples/']));
