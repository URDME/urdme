function ok = scrub_neuron(ix)
%SCRUB_NEURON Tests for NEURON Workflow.

% S. Engblom 2018-06-28

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4 @l_scrub5 @l_scrub6};
stests = {'Scrub_neuron #1' 'Scrub_neuron #2' 'Scrub_neuron #3' ...
          'Scrub_neuron #4' 'Scrub_neuron #5' 'Scrub_neuron #6'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SCRUB_NEURON (WORKFLOWS/NEURON)',ftests(ix), ...
             stests(ix));

%-------------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 ION_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
ion_run;

unix('rm mexssa*.mex*'); % cleanup
unix('rm mexuds*.mex*');
unix('rm ../../workflows/neuron/channel/lib/HHKNa_prop.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 NEURON_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
neuron_run;

unix('rm mexssa*.mex*'); % cleanup
unix('rm mexuds*.mex*');
unix('rm ../../workflows/neuron/channel/lib/HHKNa_prop.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 SYNAPTIC_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
synaptic_run;

unix('rm mexssa*.mex*'); % cleanup
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 TWONEURONS_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
twoNeurons_run;

unix('rm mexssa*.mex*'); % cleanup
unix('rm mexuds*.mex*');
unix('rm ../../workflows/neuron/channel/lib/HHKNa_prop.c');
unix('rm ../../workflows/neuron/synapse/lib/AMPA_NMDAprop.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub5
%L_SCRUB5 CIRCLENEURONS_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
circleNeurons_run;

unix('rm mexssa*.mex*'); % cleanup
unix('rm mexuds*.mex*');
unix('rm ../../workflows/neuron/channel/lib/HHKNa_prop.c');
unix('rm ../../workflows/neuron/synapse/lib/AMPA_NMDAprop.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_scrub6
%L_SCRUB6 PYRAMIDETREE_RUN example.

wd = pwd;
cd ../../workflows/neuron
ok = 1;
report = 0;
plotting_off = true;
pyramideTree_run;

unix('rm mexssa*.mex*'); % cleanup
unix('rm mexuds*.mex*');
unix('rm ../../workflows/neuron/channel/lib/HHKNa_prop.c');
unix('rm ../../workflows/neuron/synapse/lib/AMPA_NMDAprop.c');
cd(wd);

%-------------------------------------------------------------------------------
