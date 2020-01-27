function [model, geom1, std] = neuron2comsol (intree, t_vec, nodename, ...
                                              color, options)
%NEURON2COMSOL Import neuronal goemetries from TREES toolbox to Comsol
%   Needs TREES Toolbox (C) 2009  Hermann Cuntz.
%
%   Usage: Start TREES Toolbox, load geometry using load_tree() Make
%   sure Comsol LiveLink interface is loaded, then type
%   neuron2comsol()

% A. Senek 2017-07-11 (Revision)
% P. Bauer 2013-09-19
% Based on function x3d_tree from the TREES Toolbox 1.15.
% Original code by Friedrich Forstner 12 December 2008

% V 0.1 - compartments are represented as 1D-lines
% V 0.2 - import of hh-currents as time dependent current sources

% comsol intialisation
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create(nodename);
geom1 = model.geom.create('geom1',3);
model.geom('geom1').angularUnit('rad');

% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree),
    tree = trees {intree};
else
    tree = intree;
end

if (~isfield (tree, 'X')) || (~isfield (tree, 'Y'))
    [xdend tree] = xdend_tree (intree);
end

idpar = idpar_tree (intree); % vector containing index to direct parent
N =     size (idpar, 1); % number of nodes in tree

if (nargin < 4)||isempty(color),
    color = [0 0 0]; % {DEFAULT color: black}
end;
if size (color, 1) == 1,
    color = repmat (color, N, 1);
end

if (nargin < 5)||isempty(DD),
    DD = [0 0 0]; % {DEFAULT 3-tupel: no spatial displacement from the root}
end
if length(DD)<3,
    DD = [DD zeros(1, 3 - length (DD))]; % append 3-tupel with zeros
end

if (nargin < 6)||isempty(ipart),
    ipart = (1 : N)'; % {DEFAULT index: select all nodes/points}
end

if (nargin < 7)||isempty(options),
    options = '-w'; % {DEFAULT: waitbar}
end

if isfield (tree, 'D'),
    D = tree.D; % local diameter values of nodes on tree
else
    D = ones (N, 1); % if values don't exist fill with diameter = 1um
end

if strfind (options, '-thin'),
    D = ones (N, 1); % thin diameter option: all nodes 1um diameter
end
if strfind (options, '-thick'),
    D = D + 3; % thick diameter option: all nodes + 3um diameter
end

XYZ   = [tree.X tree.Y tree.Z]; % node coordinates
vXYZ  = XYZ - XYZ (idpar, :); % edge vectors
vnorm = sqrt (sum (vXYZ.^2, 2)); % norm (length) of all vectors

% raw compartments
rawComps = [zeros(N, 1) vnorm zeros(N, 1)];
% calculate compartment rotation:
% cross product to get rotation axis
rotX = (rawComps (:, 2) .* vXYZ (:, 3)) - (rawComps (:, 3) .* vXYZ (:, 2));
rotY = (rawComps (:, 3) .* vXYZ (:, 1)) - (rawComps (:, 1) .* vXYZ (:, 3));
rotZ = (rawComps (:, 1) .* vXYZ (:, 2)) - (rawComps (:, 2) .* vXYZ (:, 1));
get_zero = find (rotX == 0 & rotY == 0 & rotZ == 0);
rotX (get_zero) = 0; rotY (get_zero) = 1; rotZ (get_zero) = 0;
rotaxis = [rotX rotY rotZ];
% normalize axis
rotaxis_norm = sqrt (sum (rotaxis.^2, 2)); % norm (length) of all vectors
warning ('off', 'MATLAB:divideByZero'); rotaxis = rotaxis ./ repmat (rotaxis_norm, 1, 3);
% derive angle between rotation axis and compartment ground line
cproduct = zeros (N, 1);
for ward = 1 : N,
    cproduct (ward) = rawComps (ward, :) * vXYZ (ward, :)';
end
rotangle = acos (cproduct ./ (vnorm.^2));
warning ('on', 'MATLAB:divideByZero');
% avoid NANs
rotangle (isnan (rotangle)) = 0;
%rotationMatrix    = [rotaxis, rotangle];
%translationMatrix = repmat (DD, N, 1) + XYZ (idpar, :) + 1/2 * vXYZ;
%heightVector      = vnorm;
%radiusVector      = D/2;
%sphereMat = zeros(1,length(tree.dA));

%set geometric unity to um
geom1.lengthUnit([native2unicode(hex2dec('00b5'), 'Cp1252') 'm']);

if strfind (options, '-w'), % waitbar option: initialization
    HW = waitbar (0, 'creating geometry ...');
    set (HW, 'Name', 'writing to comsol');
end
geomCnt=1;
for ward = 2 : N,
    if strfind (options, '-w'), % waitbar option: update
        waitbar (ward ./ N, HW);
    end
    %if ipart (ward) % just to filter through function argument
    %geomCnt=ward-1;
    
    geom1.feature.create(['c',num2str(geomCnt)],'BezierPolygon');
    geom1.feature(['c',num2str(geomCnt)]).set('degree',1);
    geom1.feature(['c',num2str(geomCnt)]).set('type','open');
    geom1.feature(['c',num2str(geomCnt)]).set('createselection','on');
    
    fP=XYZ (findPredecessor(ward,tree.dA), :);
    sP=XYZ (ward,:);
    
    geom1.feature(['c',num2str(geomCnt)]).set('p', {num2str(fP(1)) num2str(sP(1)); num2str(fP(2)) num2str(sP(2)); num2str(fP(3)) num2str(sP(3))});
    geom1.feature(['c',num2str(geomCnt)]).set('w', {'1' '1'});
    
    disp(num2str(ward));
    disp(rotangle(ward, :));
    
    geomCnt=geomCnt+1;
    
    %end
end
if strfind (options, '-w'), % waitbar option: close
    close (HW);
end

if strfind (options, '->')
    if ispc,        % this even calls the file directly (only windows)
        winopen ([path tname]);
    end
end

%generate sphere around geometry
geom1.feature.create('sph1', 'Sphere');
geom1.feature('sph1').set('r', '2000');
geom1.feature('sph1').set('createselection','on');
geom1.runAll;

%create time dependent solver
model.physics.create('ec', 'ConductiveMedia', 'geom1');
model.study.create('std1');
model.study('std1').feature.create('time', 'Transient');
model.study('std1').feature('time').activate('ec', true);
std = model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'time');
model.sol('sol1').feature.create('t1', 'Time');
dt = t_vec(2) - t_vec(1);
range =  ['range(' num2str(t_vec(1)) ',' num2str(dt) ',' num2str(t_vec(end)) ')'];
model.sol('sol1').feature('t1').set('tlist', range);
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.create('i1', 'Iterative');
model.sol('sol1').feature('t1').feature('i1').set('linsolver', 'cg');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'i1');
model.sol('sol1').feature('t1').feature('i1').feature.create('mg1', 'Multigrid');
model.sol('sol1').feature('t1').feature('i1').feature('mg1').set('prefun', 'amg');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runFromTo('st1', 'v1');
model.sol('sol1').feature.create('td1', 'TimeDiscrete');
model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature('td1').set('control', 'user');
range2 =  ['range(' num2str(t_vec(1)) ',' num2str(dt)./500 ',' num2str(t_vec(end)) ')'];
model.sol('sol1').feature('td1').set('tlist', range2);
model.study('std1').feature('time').set('tunit', 'ms'); 
model.sol('sol1').feature('td1').set('control', 'time');

%   Define getCurrent as external matlab function
model.func.create('extm1', 'MATLAB');
model.func('extm1').setIndex('funcs', 'getCurrent', 0, 0);
model.func('extm1').setIndex('funcs', 'c,t', 0, 1);

%   Generate current sources for all compartments
if strfind (options, '-w'), % waitbar option: initialization
    HW = waitbar (0, 'initializing current sources ...');
    set (HW, 'Name', 'writing to comsol');
end
geomCnt=1;
for ward = 2 : N
    if strfind (options, '-w'), % waitbar option: update
        waitbar (ward ./ N, HW);
    end
    model.physics('ec').feature.create(['lcs',num2str(geomCnt)], 'LineCurrentSource', 1);
    model.physics('ec').feature(['lcs',num2str(geomCnt)]).selection.named(['geom1_c',num2str(geomCnt),'_edg']);
    model.physics('ec').feature(['lcs',num2str(geomCnt)]).set('Qjl', 1,['getCurrent(',num2str(ward),',t)']);
    geomCnt=geomCnt+1;
end
if strfind (options, '-w'), % waitbar option: close
    close (HW);
end

%%   Add model material properties ~water
model.material.create('mat1', 'Common', 'mod1');
model.material('mat1').label('Water, liquid');
model.material('mat1').set('family', 'water');
model.material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T[1/K])[Pa*s]');
model.material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.0');
model.material('mat1').propertyGroup('def').set('electricconductivity', '5.5e-6[S/m]');
model.material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T[1/K])[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('density', 'rho(T[1/K])[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', 'k(T[1/K])[W/(m*K)]');
model.material('mat1').propertyGroup('def').set('soundspeed', 'cs(T[1/K])[m/s]');
model.material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.material('mat1').propertyGroup('def').func('eta').set('funcname', 'eta');
model.material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('eta').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('eta').set('pieces', {'273.15' '413.15' '1.3799566804-0.021224019151*T^1+1.3604562827E-4*T^2-4.6454090319E-7*T^3+8.9042735735E-10*T^4-9.0790692686E-13*T^5+3.8457331488E-16*T^6'; '413.15' '553.75' '0.00401235783-2.10746715E-5*T^1+3.85772275E-8*T^2-2.39730284E-11*T^3'});
model.material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.material('mat1').propertyGroup('def').func('Cp').set('funcname', 'Cp');
model.material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('Cp').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('Cp').set('pieces', {'273.15' '553.75' '12010.1471-80.4072879*T^1+0.309866854*T^2-5.38186884E-4*T^3+3.62536437E-7*T^4'});
model.material('mat1').propertyGroup('def').func.create('rho', 'Piecewise');
model.material('mat1').propertyGroup('def').func('rho').set('funcname', 'rho');
model.material('mat1').propertyGroup('def').func('rho').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('rho').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('rho').set('pieces', {'273.15' '553.75' '838.466135+1.40050603*T^1-0.0030112376*T^2+3.71822313E-7*T^3'});
model.material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.material('mat1').propertyGroup('def').func('k').set('funcname', 'k');
model.material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.material('mat1').propertyGroup('def').func('k').set('extrap', 'constant');
model.material('mat1').propertyGroup('def').func('k').set('pieces', {'273.15' '553.75' '-0.869083936+0.00894880345*T^1-1.58366345E-5*T^2+7.97543259E-9*T^3'});
model.material('mat1').propertyGroup('def').func.create('cs', 'Interpolation');
model.material('mat1').propertyGroup('def').func('cs').set('sourcetype', 'user');
model.material('mat1').propertyGroup('def').func('cs').set('source', 'table');
model.material('mat1').propertyGroup('def').func('cs').set('funcname', 'cs');
model.material('mat1').propertyGroup('def').func('cs').set('table', {'273' '1403';  ...
'278' '1427';  ...
'283' '1447';  ...
'293' '1481';  ...
'303' '1507';  ...
'313' '1526';  ...
'323' '1541';  ...
'333' '1552';  ...
'343' '1555';  ...
'353' '1555';  ...
'363' '1550';  ...
'373' '1543'});
model.material('mat1').propertyGroup('def').func('cs').set('interp', 'piecewisecubic');
model.material('mat1').propertyGroup('def').func('cs').set('extrap', 'const');
model.material('mat1').propertyGroup('def').addInput('temperature');
model.material('mat1').set('family', 'water');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'1'});

%%  Create mesh     
mesh = model.mesh.create('mesh1', 'geom1');
mesh.autoMeshSize(4);
mesh.run;
%       % Plot the second mesh case with a colored element:
%       mphmesh(model,'mesh2','meshcolor','r');


%   Add ground for model 
model.geom('geom1').create('pt1', 'Point');
model.geom('geom1').feature('pt1').setIndex('p', '0', 0);
model.geom('geom1').feature('pt1').setIndex('p', '0', 1);
model.geom('geom1').feature('pt1').setIndex('p', '1500', 2);

model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.geom('geom1').selection('csel1').label('pnt1');
model.geom('geom1').feature('pt1').set('contributeto', 'csel1');

geom1.feature('pt1').set('createselection','on');
model.physics('ec').create('gnd1', 'Ground', 0);
model.physics('ec').feature('gnd1').selection.named('geom1_csel1_pnt');

%%  Result
model.result.dataset.create('cpl1', 'CutPlane');
model.result.dataset('cpl1').set('quickplane', 'xy');
model.result.dataset('cpl1').set('quickz', '-20');
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').feature('surf1').set('data', 'cpl1');
model.result('pg1').feature('surf1').set('colortablesym', 'on');
