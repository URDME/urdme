%STARTUP Startup for URDME.

% S. Engblom 2017-02-14

% link = location of this startup.m
link = mfilename('fullpath');
link = link(1:end-7); % s-t-a-r-t-u-p is 7 chars

% path to urdme/ folders
addpath(genpath([link 'urdme/']));
addpath(genpath([link 'workflows/']));