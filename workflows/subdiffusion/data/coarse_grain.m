%COARSE_GRAIN Coarse-grained internal states model of crowding.

% S. Engblom 2017-05-30

% Statistics for spherical crowders
clc
clear all
close all

 
R = 1e-1; % Radius crowder
r = 1e-1; % Radius moving molecule

 
repeats = 1e2;

 
% expected first exit times
E = zeros(1,repeats);
LAMBDA = zeros(repeats,2);
crowdcount = 0;

 
for Crowding = 0.1:0.1:0.4
    Crowding
    NCrowders = ceil(Crowding/(R^2));
    crowdcount = crowdcount+1;
fprintf(1,'Repeat: ');
for count=1:repeats

    
    % display counter
    if count>1
      for j=0:log10(count-1)
          fprintf('\b'); % delete previous counter display
      end
    end
    fprintf(1,'%d',count);

    

    
    Positions = zeros(NCrowders,2);
    for k = 1:NCrowders
        accept = 0;
        while ~accept
            x = 2*rand-1;
            y = 2*rand-1;
            if k == 1
                if (sqrt(x^2+y^2)<(1+R+r))  && (sqrt(x^2+y^2)>(R+r))
                    accept = 1;
                    Positions(k,:) = [x,y];
                end
            else
                dist = sqrt((Positions(1:k-1,1)-x).^2+(Positions(1:k-1,2)-y).^2);
                if  (sqrt(x^2+y^2)<(1+R+r)) && (sum(dist >= 2*R)==k-1)  && (sqrt(x^2+y^2)>(R+r))
                    accept = 1;
                    Positions(k,:) = [x,y];
                end
            end
        end
    end

    
    import com.comsol.model.*
    import com.comsol.model.util.*

 
    model = ModelUtil.create('Model');

 
    model.modelPath('/Users/linameinecke/Documents/phd/ComsolPhysics/Crowding');

 
    model.comments(['untitled\n\n']);

 
    model.modelNode.create('comp1');

 
    model.geom.create('geom1', 2);

 
    model.mesh.create('mesh1', 'geom1');

 
    model.physics.create('poeq', 'PoissonEquation', 'geom1');

 
    model.study.create('std1');
    model.study('std1').create('stat', 'Stationary');
    model.study('std1').feature('stat').activate('poeq', true);

 
    % Create Exit circle
    model.geom('geom1').run('');
    model.geom('geom1').feature.create('c1', 'Circle');
    model.geom('geom1').feature('c1').set('type', 'solid');
    model.geom('geom1').feature('c1').set('base', 'center');
    model.geom('geom1').feature('c1').set('pos', {'0' '0'});
    model.geom('geom1').feature('c1').set('r', '1');
    model.geom('geom1').run('c1');
    model.geom('geom1').run('c1');
    % Create point at origin
    model.geom('geom1').feature.create('pt1', 'Point');
    model.geom('geom1').feature('pt1').set('p', {'0' '0'});
    model.geom('geom1').run('pt1');
    model.geom('geom1').run;

 
    % Mesh size
    model.mesh('mesh1').autoMeshSize(6);
    % order of FEM basis functions
    model.physics('poeq').prop('ShapeProperty').set('order', '1');

 
    % Select segments for homogenous Dirchlet b.c.
    model.geom('geom1').create('sel1', 'ExplicitSelection');
    model.geom('geom1').feature('sel1').selection('selection').init(1);
    model.geom('geom1').feature('sel1').selection('selection').set('c1', [1 2 3 4]);
    % lower left segment -> left hop
    model.geom('geom1').create('sel2', 'ExplicitSelection');
    model.geom('geom1').feature('sel2').selection('selection').init(1);
    model.geom('geom1').feature('sel2').selection('selection').set('c1', [1]);
    % lower right segment -> down hop
    model.geom('geom1').create('sel3', 'ExplicitSelection');
    model.geom('geom1').feature('sel3').selection('selection').init(1);
    model.geom('geom1').feature('sel3').selection('selection').set('c1', [2]);
    % upper right segment -> right hop
    model.geom('geom1').create('sel4', 'ExplicitSelection');
    model.geom('geom1').feature('sel4').selection('selection').init(1);
    model.geom('geom1').feature('sel4').selection('selection').set('c1', [3]);
    % upper left segment -> up hop
    model.geom('geom1').create('sel5', 'ExplicitSelection');
    model.geom('geom1').feature('sel5').selection('selection').init(1);
    model.geom('geom1').feature('sel5').selection('selection').set('c1', [4]);
    % select point
    model.geom('geom1').create('sel6', 'ExplicitSelection');
    model.geom('geom1').feature('sel6').selection('selection').init(0);
    model.geom('geom1').feature('sel6').selection('selection').set('pt1', [1]);

 
    % Set homogeneous Dirichlet on all boundaries
    model.physics('poeq').feature.create('dir1', 'DirichletBoundary', 1);
    model.physics('poeq').feature('dir1').selection.named('geom1_sel1');

 
    model.geom('geom1').run('sel1');
    model.geom('geom1').run('sel2');
    model.geom('geom1').run('sel3');
    model.geom('geom1').run('sel4');
    model.geom('geom1').run('sel5');
    model.geom('geom1').run('sel6');

 
    Names = [];

 
    if NCrowders >0
        % Create crowding molecule
        for k = 1:NCrowders

 
            x = Positions(k,1);
            y = Positions(k,2);

 
            CircleName = sprintf('c%d',k+1);
            Names = [Names, CircleName, ' '];

 
            % Create crowding molecule
            model.geom('geom1').feature.create(CircleName, 'Circle');
            model.geom('geom1').feature(CircleName).set('type', 'solid');
            model.geom('geom1').feature(CircleName).set('base', 'center');
            model.geom('geom1').feature(CircleName).set('pos', {num2str(x) num2str(y)});
            model.geom('geom1').feature(CircleName).set('r', num2str(R+r));
            model.geom('geom1').run(CircleName);
            model.geom('geom1').run(CircleName);
        end

 
        Names = strsplit(Names(1:end-1));
        % Create geometry = Exit circle - crowders
        model.geom('geom1').create('dif1', 'Difference');
        model.geom('geom1').feature('dif1').selection('input').set({'c1'});
        % Cut out crowders
        model.geom('geom1').feature('dif1').selection('input2').set(Names);
    end

 
%         figure
%         % plot geometry with/without boundary lables
%         mphgeom(model,'geom1','Edgelabels','off','Vertexlabels','off','facelabels','off')
% %         % plot molecule radius
% %         viscircles(PositionsLocal,(RCrowder)*ones(NCrowdersLocal,1));
%         % write boundary segments
%         model.selection
%         mphgetselection(model.selection('geom1_sel1'))
%         mphgetselection(model.selection('geom1_sel2'))
%         mphgetselection(model.selection('geom1_sel3'))
%         mphgetselection(model.selection('geom1_sel4'))
%         mphgetselection(model.selection('geom1_sel5'))
%         mphgetselection(model.selection('geom1_sel6'))

 
    model.geom('geom1').run;

    
    % Select domain without excluded volume
    model.geom('geom1').create('ballsel1', 'BallSelection');
    model.geom('geom1').feature('ballsel1').set('posx', '0');
    model.geom('geom1').feature('ballsel1').set('posy', '0');
    model.geom('geom1').feature('ballsel1').set('r', num2str((R+r)/1000));
    model.geom('geom1').runPre('fin');
    model.geom('geom1').run;
    model.physics('poeq').selection.named('geom1_ballsel1');

    
    center_domain = model.selection('geom1_ballsel1').entities(2); % number central subdomain

    
    upDown = model.geom('geom1').getUpDown;
    domain = cell(1,max(upDown(:))+1);
    for i = 1:length(upDown)
    domain{upDown(1,i)+1}=[domain{upDown(1,i)+1},i];
    domain{upDown(2,i)+1}=[domain{upDown(2,i)+1},i];
    end

     
    if isempty(center_domain) % node cut out
            E(count) = inf;
            disp('Node cut out')

    
    else
        center_boundaries = domain{center_domain+1}; % boundaries connected to central subdomain
        dirichlet_boundaries = model.selection('geom1_sel1').entities(1);

    
        if isempty(intersect(center_boundaries,dirichlet_boundaries)) %central subdomain not connected to dirichlet bc
            E(count) = inf;

        
        else
            model.sol.create('sol1');
            model.sol('sol1').study('std1');

 
            model.study('std1').feature('stat').set('notlistsolnum', 1);
            model.study('std1').feature('stat').set('notsolnum', '1');
            model.study('std1').feature('stat').set('listsolnum', 1);
            model.study('std1').feature('stat').set('solnum', '1');

 
            model.sol('sol1').create('st1', 'StudyStep');
            model.sol('sol1').feature('st1').set('study', 'std1');
            model.sol('sol1').feature('st1').set('studystep', 'stat');
            model.sol('sol1').create('v1', 'Variables');
            model.sol('sol1').feature('v1').set('control', 'stat');
            model.sol('sol1').create('s1', 'Stationary');
            model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
            model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
            model.sol('sol1').feature('s1').feature.remove('fcDef');
            model.sol('sol1').attach('std1');

 
            model.result.create('pg1', 2);
            model.result('pg1').set('data', 'dset1');
            model.result('pg1').create('surf1', 'Surface');
            model.result('pg1').feature('surf1').set('expr', 'u');

 
            model.sol('sol1').runAll;

 
            model.result('pg1').run;

 
            % Evalute E at point p1
            model.result.numerical.create('pev1', 'EvalPoint');
            model.result.numerical('pev1').selection.named('geom1_sel6');
            model.result.table.create('tbl1', 'Table');
            model.result.table('tbl1').comments('Point Evaluation 1 (u)');
            model.result.numerical('pev1').set('table', 'tbl1');
            model.result.numerical('pev1').setResult;
            % Expected exit time
            e = model.result.table('tbl1').getReal();

 
            E(count) = e*4; % E/E_0 (E_0 = 1/4)
        end
    end
end

 
lambda = 1./E;
LAMBDA(:,crowdcount) = lambda;

 
end
out = model;
