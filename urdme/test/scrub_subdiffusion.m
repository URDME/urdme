function ok = scrub_subdiffusion(ix)
%SCRUB_SUBDIFFUSION Tests for SUBDIFFUSION workflow.

% S. Engblom 2024-05-06

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4};
stests = {'Scrub_subdiffusion #1' 'Scrub_subdiffusion #2' ...
          'Scrub_subdiffusion #3' 'Scrub_subdiffusion #4'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SCRUB_SUBDIFFUSION (WORKFLOWS/SUBDIFFUSION)', ...
             ftests(ix),stests(ix));

%-------------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 Test of pure1D_run.

ok = 1;

wd = pwd;
cd ../../workflows/subdiffusion
ok = 1;
plotting_off = true;
pure1D_run;
unix('rm mexnsm*.mex*'); % cleanup
unix('rm mexuds*.mex*');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 Test of pure2D_run.

ok = 1;

wd = pwd;
cd ../../workflows/subdiffusion
ok = 1;
report = 0;
plotting_off = true;
pure2D_run;
unix('rm mexnsm*.mex*'); % cleanup
unix('rm mexuds*.mex*');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 Test of bimolecular2D_run.

ok = 1;

wd = pwd;
cd ../../workflows/subdiffusion
ok = 1;
report = 0;
plotting_off = true;
t = linspace(0,2,11);

Csave = cell(1,5);
for Kcase = 1:5
  bimolecular2D_run
  Csave{Kcase} = C;
end

unix('rm mexnsm*.mex*'); % cleanup
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 Test of mindcde3D_run.

ok = 1;

wd = pwd;
cd ../../workflows/subdiffusion
ok = 1;
report = 0;
tspan = 0:20;
mincde3D_run;
unix('rm mexnsm*.mex*'); % cleanup
unix('rm mincde_subdiff.c'); 
cd(wd);
%-------------------------------------------------------------------------------


