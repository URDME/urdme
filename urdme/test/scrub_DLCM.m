function ok = scrub_DLCM(ix)
%SCRUB_DLCM Tests for DLCM Workflow.

% E. Blom 2024-05-20 (solver scrubs)
% S. Engblom 2023-09-19

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4 @l_scrub5 @l_scrub6 ...
          @l_scrub7 @l_scrub8};
stests = {'Scrub_DLCM #1' 'Scrub_DLCM #2' 'Scrub_DLCM #3' ...
          'Scrub_DLCM #4' 'Scrub_DLCM #5' 'Scrub_DLCM #6' ...
          'Scrub_DLCM #7' 'Scrub_DLCM #8'};
if nargin == 0, ix = 1:size(ftests,2); end
str = 'SCRUB_DLCM';
ok = runtest(str,ftests(ix), stests(ix));

%-------------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 Relaxation from initially crowded population of cells.
% tests basic pressure mechanics and that solver runs without user defining
% surface tension, curvature operators, and internal reaction mechanics

ok = 1;
wd = pwd;
cd ../../workflows/DLCM/experiments/1_basic_test
ok = 1;
Tend = 5000;
basic_test

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm basic_test_outer.c');
cd(wd);

% check correctness
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...    % get U history
size(umod.U,2),size(umod.U,3));
ok = ok && (sum(U(:,end)==1)>0.85*sum(U(:,end)>0)); % final size
ok = ok && ~any(sum(U) ~= 544);                 % conservation of cells

%-------------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 Basic_test in 3D with surface tension

ok = 1;
if exist('fegeometry.m','file') % 3D PDE Toolbox required
  wd = pwd;
  cd ../../workflows/DLCM/experiments/1_basic_test/
  Tend = 1000;
  basic_test_3D

  unix('rm mexdlcm*.mex*'); % cleanup
  unix('rm basic_test_3D_outer.c');
  cd(wd);

  U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
  size(umod.U,2),size(umod.U,3));
  ok = ok && ~any(sum(U) ~= 1428);       % mass conservation
  Adof = find(U(:,end));
  ok = ok && (abs(max(P(1,Adof).^2 + P(2,Adof).^2) - 0.5430) < 0.05);
end

%-------------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 Avascular tumor model
% tests surface tension and reaction events (phenotype switching and birth)

ok = 1;
wd = pwd;
cd ../../workflows/DLCM/experiments/2_avascular_tumour
ok = 1;
Tend = 2000;
avascular_tumor

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm avascular_tumor_outer.c');
cd(wd);

% check correctness
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
Adof = find(U(:,end));
X = umod.U(:,:,:);
ok = ok && (abs(max(P(1,Adof).^2 + P(2,Adof).^2) - 0.0469) < 0.01);

ok = ok && (abs(sum(X(1,:,end))-3) < 1);         % nr degrading cells
ok = ok && (max(abs(X(3,:,end))) < 0.005);        % pressure levels
ok = ok && (abs(min(X(4,:,end))-0.9461) < 0.005);  % oxygen levels

%-------------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 Chemotaxis model using Darcy's law _and_ chemotaxis upwards a
% chemical gradient

ok = 1;
wd = pwd;
cd ../../workflows/DLCM/experiments/3_gradient_growth/
ok = 1;
Tend = 100;
chemotaxis

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm chemotaxis_outer.c');
cd(wd);

% check correctness
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
Adof = find(U(:,end));
ok = ok && (abs(min(P(1,Adof)) + 0.255) < 0.02);
ok = ok && ~any(diff(sum(U)));          % conservation of number of cells

%-------------------------------------------------------------------------------
function ok = l_scrub5
%L_SCRUB5 Chemotaxis model in 3D. Testing several sources of
% movement and 3D space.

ok = 1;
if exist('fegeometry.m','file') % 3D PDE Toolbox required
  wd = pwd;
  cd ../../workflows/DLCM/experiments/3_gradient_growth
  Tend = 100;
  chemotaxis3D

  unix('rm mexdlcm*.mex*'); % cleanup
  unix('rm chemotaxis3D_outer.c');
  cd(wd);

  % check correctness
  U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
  size(umod.U,2),size(umod.U,3));
  Adof = find(U(:,end));
  ok = ok && ~any(sum(U) ~= 711);       % mass conservation
  ok = ok && abs(var(P(1,Adof).^2 + P(2,Adof).^2 + P(3,Adof).^2)-0.0018) ...
    < 0.0002; % spread of cells
end

%-------------------------------------------------------------------------------
function ok = l_scrub6
%L_SCRUB6 continuous Delta-Notch dynamics using Euler Forward

wd = pwd;
cd ../../workflows/DLCM/experiments/4_delta_notch
ok = 1;
Tend = 100;
internal_state = 'cont';  % impose continuous internal states
seed = 123;
delta_notch

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm mexnsm*.mex*');
unix('rm delta_notch_outer.c');
unix('rm delta_notch_inner.c');
cd(wd);

% check correctness
idx = find(strcmp([umod.solverargs(:)], 'mumod')) + 1; % find mumod
data = reshape(umod.solverargs{idx}.U(1,:,:), ...
size(umod.U,2), size(umod.U,3));             % = mumod.u0(1,:,:)
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
Adof = find(U(:,end));
% delta mean within certain value
ok = ok && abs(mean(data(Adof,end)) - 29) < std(data(Adof,end));

%-------------------------------------------------------------------------------
function ok = l_scrub7
%L_SCRUB7 Discrete Delta-Notch signaling using SSA solver inside DLCM.

wd = pwd;
cd ../../workflows/DLCM/experiments/4_delta_notch
ok = 1;
Tend = 100;
seed = 123;
delta_notch               % discrete internal states are default

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm mexnsm*.mex*');
unix('rm mexuds*.mex*');
unix('rm delta_notch_outer.c');
unix('rm delta_notch_inner.c');
cd(wd);

% check correctness
idx = find(strcmp([umod.solverargs(:)],'mumod')) + 1; % find mumod
data = reshape(umod.solverargs{idx}.U(1,:,:), ...
size(umod.U,2), size(umod.U,3));             % = mumod.u0(1,:,:)
U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
Adof = find(U(:,end));
% delta mean within certain value
ok = ok && abs(mean(data(Adof,end)) - 29) < std(data(Adof,end));

%-------------------------------------------------------------------------------
function ok = l_scrub8
%L_SCRUB8 Cell sorting by Y-L pressure drop; 2 cell types.

wd = pwd;
cd ../../workflows/DLCM/experiments/9_cell_sorting
ok = 1;
Tend = 100;
warning('off','rparse:ghost_species') % cell types with no reaction events
cell_sorting
warning('on','rparse:ghost_species')

unix('rm mexdlcm*.mex*'); % cleanup
unix('rm cell_sorting_outer.c');
cd(wd);

U = reshape(sum(umod.U(1:ntypes,:,:),1), ...
size(umod.U,2),size(umod.U,3));
Adof = find(U(:,end));
ok = ok && ~any(sum(U) ~= 485);       % mass conservation
% cell type conservation
ok = ok && ~sum(sum(umod.U(1:ntypes,:,end),2) - ...
  sum(umod.U(1:ntypes,:,1),2));

%-------------------------------------------------------------------------------
