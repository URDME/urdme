function animate_growth(filename, modeltype, save_animation)
% create a GIF animation from file


loadsim = ...%["DLCM_example2.mat";
           [filename];
       
%modeltype = ["PDE"] or ["DLCM"];
       
nsims = numel(loadsim);
Usave_all = cell(nsims, 1);
Oxysave_all = cell(nsims, 1);
kappa_prol_all = cell(nsims, 1);
kappa_death_all = cell(nsims, 1);
R_all = cell(nsims, 1);
V_all = cell(nsims, 1);
 
% load animation data
for n = 1:nsims
    load(loadsim(n), 'Usave','Oxysave','all_parameters','R','V')
    Usave_all{n} = Usave; Oxysave_all{n} = Oxysave; % cell and oxygen data
    kappa_prol_all{n} = all_parameters(2);      % oxygen parameters
    kappa_death_all{n} = all_parameters(4);
    R_all{n} = R;                               % geometry
    V_all{n} = V;
    timesteps_all(n) = numel(Usave);
end

% stop if not all animations are equal length
if std(timesteps_all) ~= 0
    error('animations are not the same length.')
end

timesteps = timesteps_all(1);

VUsave_all = cell(nsims, 1);
% get visited voxels
for n = 1:nsims
    VUsave_all{n} = cell(1,timesteps);
    VU_tmp = (Usave_all{n}{1} ~= 0);
    VUsave_all{n}{1} = VU_tmp;
    for i = 2:timesteps
        VU_tmp = VU_tmp | Usave_all{n}{i};
        VUsave_all{n}{i} = VU_tmp;
    end
end

% population appearance
tpf = 5;    % time steps per frame
im = 1; % individual index for M struct to make it proper size (if tpf>1)

M = struct('cdata',{},'colormap',{});
figure(3), clf,
set(gcf,'Position',[100 100 480 480]);
for i = 2:tpf:timesteps
  
  for n = 1:nsims % for all simulations
      
      subplot(1, 1, n)
      
      if modeltype(n) == "DLCM"
      pdof = find(Usave_all{n}{i} > 0 & (Oxysave_all{n}{i} > kappa_prol_all{n})); % prolif. region
      qdof = find(Usave_all{n}{i} > 0 & (Oxysave_all{n}{i} <= kappa_prol_all{n}) ...
                               & (Oxysave_all{n}{i} > kappa_death_all{n})); % quiesc. region
      ndof = find(Usave_all{n}{i} == - 1);  % region of necrotic cells (note difference with above)
      ddof = find(Usave_all{n}{i} > 0 & (Oxysave_all{n}{i} <= kappa_death_all{n})); % living, but dying
      patch('Faces',R_all{n},'Vertices',V_all{n},'FaceColor',[1 1 1], ...
            'EdgeColor','none');
      
      % all visited voxels
      outdof = find(VUsave_all{n}{i} & Usave_all{n}{i} == 0);
        
      % light grey background
      patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
      hold on,
      axis(0.75.*[-1 1 -1 1]); axis square, axis off
       
      c = [0.6 0.6 0.6];   % darker grey than background
      patch('Faces',R(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');
      
      %ii = find(Usave{i} == 1);
      single = patch('Faces',R_all{n}(qdof,:),'Vertices',V_all{n}, ...
            'FaceColor',graphics_color('bluish green'), 'EdgeColor','none');

      %ii = find(Usave{i} == 2);
      double = patch('Faces',R_all{n}(pdof,:),'Vertices',V_all{n}, ...
            'FaceColor',graphics_color('vermillion'), 'EdgeColor','none');

      %ii = find(Usave{i} == -1);
      dead = patch('Faces',R_all{n}(ndof,:),'Vertices',V_all{n}, ...
            'FaceColor',[0 0 0], 'EdgeColor','none');

      dying = patch('Faces',R_all{n}(ddof,:),'Vertices',V_all{n}, ...
            'FaceColor',[0.5 0.5 0.5], 'EdgeColor','none'); 
 
      % mark double
      ii = find(Usave_all{n}{i} == 2);
      patch('Faces',R(ii,:),'Vertices',V, ...
      'FaceColor',[0 0 0], 'FaceAlpha', 0.45, 'EdgeColor','none');
      end
     
     % PDE
     if modeltype(n) == "PDE"
     cutoff_bdof = all_parameters(6);
     pdof = find(Usave_all{n}{i} > cutoff_bdof & ... % non-boundary prolif. dofs
     (Oxysave_all{n}{i} > kappa_prol_all{n}));  
    qdof = find(Usave_all{n}{i} > cutoff_bdof & ... % quiescent dofs
     (Oxysave_all{n}{i} <= kappa_prol_all{n}) & (Oxysave_all{n}{i} >= kappa_death_all{n}));
    ndof = find(Usave_all{n}{i} > cutoff_bdof & ... % 'necrotic' dof
     (Oxysave_all{n}{i} < kappa_death_all{n})); 

    % nonzero density outside tumor
    outdof = find(Usave_all{n}{i} & Usave_all{n}{i} <= cutoff_bdof); 
    
    % light grey background
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis(0.75.*[-1 1 -1 1]); axis square, axis off
    
    % plot regions
    c = graphics_color('vermillion');
    patch('Faces',R_all{n}(pdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none');     
    c = graphics_color('bluish green');
    patch('Faces',R_all{n}(qdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0 0 0];
    patch('Faces',R_all{n}(ndof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
    c = [0.6 0.6 0.6];   % darker grey than background
    patch('Faces',R_all{n}(outdof,:),'Vertices',V,'FaceVertexCData',c, ... 
     'FaceColor','flat', 'EdgeColor','none'); 
     end
 end
  
  %subplot(2, 2, 1)
  %legend([single, double, dead],'quiescent','proliferating', 'necrotic');
  
  %title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
  drawnow;

  %set(gcf, 'Position', [100,100, 300, 300])
  M(im) = getframe(gcf);
  %writeVideo(v,M(im))
  im = im+1;
end

% save animaiton as gif
if save_animation
    movie2gif(M,{M([1:2 end]).cdata},'PDE_exp15_5.gif', ...
          'delaytime',0.1,'loopcount',0);
end
end