function ok = examples(ix)
%EXAMPLES Tests for URDME.

% S. Engblom 2017-05-09 (removed Comsol dependencies)
% S. Engblom 2017-03-17 (morphogenesis3D)
% S. Engblom 2017-03-14 (mincde)
% S. Engblom 2017-03-09 (morphogenesis)
% S. Engblom 2017-02-19

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_example1 @l_example2 @l_example3 @l_example4 @l_example5 ...
          @l_example6 @l_example7};
stests = {'Example #1/annihilation' 'Example #2/annihilation2D' ...
          'Example #3/basic' 'Example #4/morphogenesis2D' ...
          'Example #5/morphogenesis3D' 'Example #6/mincde' ...
          'Example #7/Lotka-Volterra'};

if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('EXAMPLES (URDME)', ftests(ix), stests(ix));

%-------------------------------------------------------------------------------
function ok = l_example1
%L_EXAMPLE1 Annihilation.

wd = pwd;
cd ../../examples/annihilation
ok = 1;
plotting_off = true;
annihilation_run;

ma = 1e2*[5.684618722554892 5.120565229540920 0 0]';
mb = 1e2*[0 0 4.260109840319362 4.826490738522955]';
ok = ok && norm(mean_A([1 2 end-1 end])-ma,inf) < 1e-5;
ok = ok && norm(mean_B([1 2 end-1 end])-mb,inf) < 1e-5;
unix('rm mexnsm*.mex*'); % cleanup
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example2
%L_EXAMPLE2 Annihilation 2D, requires PDE Toolbox.

wd = pwd;
cd ../../examples/annihilation
ok = 1;
plotting_off = true;
annihilation2D_run;

% it is difficult to test in a good way across different platforms:
ok = ok && sum(umod.U(1:2:end,end).*(umod.sd == 3)) <= 10;
ok = ok && sum(umod.U(2:2:end,end).*(umod.sd == 1)) <= 5;
unix('rm mexnsm*.mex*'); % cleanup
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example3
%L_EXAMPLE3 Basic. This test uses a Comsol Livelink, if present.

wd = pwd;
cd ../../examples/basic
ok = 1;
plotting_off = true;
basic_run;

Z_series_fastDiff = sum(Z_series_fastDiff);
Z_series_fastDiff = Z_series_fastDiff([1:2 end-1:end]);
Z_series_slowDiff = sum(Z_series_slowDiff);
Z_series_slowDiff = Z_series_slowDiff([1:2 end-1:end]);
% it is difficult to test in a good way across different platforms:
%ok = ok && all(Z_series_fastDiff == [0    14    38    34]);
%ok = ok && all(Z_series_slowDiff == [0    11    33    34]);
% cleanup
unix('rm mexnsm*.mex*');
unix('rm basic.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example4
%L_EXAMPLE4 Morphogenesis, requires PDE Toolbox.

wd = pwd;
cd ../../examples/morphogenesis
ok = 1;
plotting_off = true;
tspan = 0:0.1:1;
morphogenesis2D_run;
% cleanup
unix('rm mexnsm*.mex*');
unix('rm schnakenberg.c');
unix('rm brusselator.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example5
%L_EXAMPLE5 Morphogenesis 3D. This test uses a Comsol Livelink, if present.

wd = pwd;
cd ../../examples/morphogenesis
ok = 1;
plotting_off = true;
tspan = 0:0.1:1;
morphogenesis3D_run;
% cleanup
unix('rm mexnsm*.mex*');
unix('rm schnakenberg.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example6
%L_EXAMPLE6 MinCDE. This test uses a Comsol Livelink, if present.

wd = pwd;
cd ../../examples/mincde
ok = 1;
plotting_off = true;
tspan = 0:0.1:1;
mincde_run;
% cleanup
unix('rm mexnsm*.mex*');
unix('rm fange.c');
cd(wd);

%-------------------------------------------------------------------------------
function ok = l_example7
%L_EXAMPLE7 Lotka-Volterra.

wd = pwd;
cd ../../examples/Lotka_Volterra
ok = 1;
plotting_off = true;
tspan = 0:0.1:10;
LV_run;
% cleanup
unix('rm mexnsm*.mex*');
cd(wd);

%-------------------------------------------------------------------------------
