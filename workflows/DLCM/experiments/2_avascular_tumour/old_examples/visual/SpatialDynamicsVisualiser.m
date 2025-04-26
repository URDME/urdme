% Plots figs and gifs of PDE tumor spatial solution
%
% Loads simulation 'PDE_exp#.mat' or 'DLCM_exp#.mat' (user input #)
% and plots solution at times `tstamps`


% E. Blom 2023-04-01

clear;

% Animate PDE or DLCM
PDE = true;
savefigs = false;

abc = ['a', 'b', 'c']; % filename ending
testnr = 2; % 1 (k=2), 2 (k=1), 15 (k=3), 16 (k=1), 20 (k=2), 22 (k=inf).

if PDE
    testdata = 'PDE_exp' + string(testnr);
    % time stamps to plot
    %tstamps = [1 101 121];
    %tstamps = [2 101 151];
    %tstamps = [1 101 801];
    tstamps = [1 94 112]; % scaled by 2.7

else
    testdata = 'DLCM_exp' + string(testnr);
    % time stamps to plot
    tstamps = [1 85 101];  % default tstamp values
    %tstamps = [1 63 76];  % if tspan -> 40
    %tstamps = [2 71 97];
end

load(testdata)

i = 1;
if PDE
% extract parameters if loading simulation
disp('Extracting parameters from loaded .mat-file')
mu_prol = all_parameters(3);
mu_death = all_parameters(5);
lambda = all_parameters(1);
kappa_prol = all_parameters(2);
kappa_death = all_parameters(4);
cutoff_bdof = all_parameters(6);
sigma = all_parameters(7);

%% visualize PDE
% get VU for each tstamp if nonexistent
if  ~exist('VUsave', 'var')
    VUsave = cell(1,numel(tspan));
    VU_tmp = (Usave{1} ~= 0);
    for i = 2:numel(tspan)
        VU_tmp = VU_tmp | Usave{i}>cutoff_bdof;
        VUsave{i} = VU_tmp;
    end
end

i=1;
for iii = tstamps
    fig=figure(i);
    set(gcf,'Position',[100 100 340/3 240/2]);
    clf,

    Umat=full(cell2mat(Usave));

    % background
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
     'EdgeColor','none');

    hold on,
    axis(0.75.*[-1 1 -1 1]); axis square, axis off
     
    if iii ~= 1
    % colour tumor volumes after growth/death/quiescence

    pdof = find(Usave{iii} > cutoff_bdof & ... % non-boundary prolif. dofs
     (Oxysave{iii} > kappa_prol));  
    qdof = find(Usave{iii} > cutoff_bdof & ... % quiescent dofs
     (Oxysave{iii} <= kappa_prol) & (Oxysave{iii} >= kappa_death));
    ndof = find(Usave{iii} > cutoff_bdof & ... % 'necrotic' dof
     (Oxysave{iii} < kappa_death)); 

    % nonzero density outside tumor
    %outdof = find(Usave{iii} & Usave{iii} <= cutoff_bdof); 
    outdof = find(VUsave{iii} & Usave{iii} <= cutoff_bdof); 

    % plot regions
    c = graphics_color('vermillion');
    patch('Faces',R(pdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');     
    c = graphics_color('bluish green');
    patch('Faces',R(qdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0 0 0];
    patch('Faces',R(ndof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0.6 0.6 0.6];   % darker grey than background
    patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 

    else
        pdof = find(Usave{iii} > cutoff_bdof);
        c = graphics_color('vermillion');
        patch('Faces',R(pdof,:),'Vertices',V,'FaceVertexCData',c, ... 
            'FaceColor','flat', 'EdgeColor','none'); 
    end

    % mark domain origo by crosshairs
    plot([0 0], [-1 1],'w')
    plot([-1 1], [0 0] ,'w')
    drawnow;
   
    % % uncomment to save the GIF
    %Tumour(iii) = getframe(gcf);
    %movie2gif(Tumour,{Tumour.cdata},'allidof_test1.gif', ...
    %       'delaytime',0.1,'loopcount',0);

    set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
    figname = 'PDE_spatialdynamics_' + ...
                  string(testnr) + string(abc(i)) + '.pdf';
    if savefigs
    	exportgraphics(gca,figname)
    end
    
    i = i+1;
end
else
    
disp('Extracting parameters from loaded .mat-file')
mu_prol = all_parameters(3);
mu_death = all_parameters(5);
lambda = all_parameters(1);
kappa_prol = all_parameters(2);
kappa_death = all_parameters(4);
mu_degrade = all_parameters(6);
sigma = all_parameters(7);

%% visualise DLCM
% get VU for each tstamp if nonexistent
if  ~exist('VUsave', 'var')
    VUsave = cell(1,numel(tspan));
    VU_tmp = (Usave{1} ~= 0);
    for i = 2:numel(tspan)
        VU_tmp = VU_tmp | Usave{i};
        VUsave{i} = VU_tmp;
    end
end

i = 1;
for iii = tstamps
    fig=figure(i);
    set(gcf,'Position',[100 100 340/3 240/2]);
    clf,

    
  %TEMP: evaluate oxygen for runs without Oxysave
  %adof = find(Usave{i}); % all filled voxels
  %sdof = find(Usave{i} > 1); % voxels with 2 cells
    % RHS source term proportional to the over-occupancy and BCs
  %Oxy = full(fsparse([extdof; adof],1, ...
                     %[ones(size(extdof)); ... 
                     % -cons*full(max(Usave{i}(adof),0).*dM(adof))], ...
                     %[size(OLa.X,1) 1]));
  %Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));
  %Oxysave{i} = Oxy;
  
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
  hold on,
  axis(0.75.*[-1 1 -1 1]); axis square, axis off

  if iii ~= 1
  pdof = find(Usave{iii} > 0 & (Oxysave{iii} > kappa_prol)); % prolif. region
  qdof = find(Usave{iii} > 0 & (Oxysave{iii} <= kappa_prol) ...
                           & (Oxysave{iii} > kappa_death)); % quiesc. region
  ndof = find(Usave{iii} == - 1);  % region of necrotic cells (note difference with above)
  ddof = find(Usave{iii} > 0 & (Oxysave{iii} <= kappa_death)); % living, but dying
  
  outdof = find(VUsave{iii} & Usave{iii} == 0);
  
  c = [0.6 0.6 0.6];   % darker grey than background
  patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');
 
  %ii = find(Usave{i} == 1);
  single = patch('Faces',R(qdof,:),'Vertices',V, ...
        'FaceColor',graphics_color('bluish green'), 'EdgeColor','none');
  
  %ii = find(Usave{i} == 2);
  double = patch('Faces',R(pdof,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'), 'EdgeColor','none');
  
  %ii = find(Usave{i} == -1);
  dead = patch('Faces',R(ndof,:),'Vertices',V, ...
        'FaceColor',[0 0 0], 'EdgeColor','none');
    
  dying = patch('Faces',R(ddof,:),'Vertices',V, ...
        'FaceColor',[0.25 0.25 0.25], 'EdgeColor','none');  
    
  % mark double
  ii = find(Usave{iii} == 2);
  patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0], 'FaceAlpha', 0.45, 'EdgeColor','none');
    
  else
      % assuming initial state is U=1 with fully oxygen saturation
        pdof = find(Usave{iii} > 0);
        double = patch('Faces',R(pdof,:),'Vertices',V, ...
        'FaceColor',graphics_color('vermillion'), 'EdgeColor','none');
  end

  % mark domain origo by crosshairs
  plot([0 0], [-1 1],'w')
  plot([-1 1], [0 0] ,'w')
    
  %legend([single, double, dead],'quiescent','proliferating', 'necrotic');
 
  drawnow;
  
  set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
  figname = 'DLCM_spatialdynamics_' + ...
                  string(testnr) + string(abc(i)) + '.pdf';
  if savefigs
    exportgraphics(gca,figname)
  end
  
  i = i+1;
end
end
 
% uncomment to save
%exportgraphics(gca,'spatialdynamics_test9_a.pdf')