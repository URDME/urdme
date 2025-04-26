%VISUALISE_AVASCULAR_TUMOR
%   Plots 5 snapshots of tumor growth from simulation data (assumes
%   certain oxygen thresholds)

% E. Blom 2025-04-17

%load tumor1.mat        % sigma > 0; tumor remains circular
load tumor2.mat         % sigma = 0; tumor breaks
U = umod.U(1:2,:,:);    % nr of cell types per voxel
Q = umod.U(3:4,:,:);    % pressure & nutrients

kappa_prol = 0.955;    % oxygen threshold for proliferation
tspan = 1:5; % times for each snapshot
for i = tspan
  figure(),
  clf,
  hold on
  patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
  'EdgeColor','none');
  Utot = sum(U(:,:,i),1);     % nr cells per voxel
  patch('Faces',R(Utot>0&Q(2,:,i)<=kappa_prol,:),'Vertices',V, ...
  'FaceColor',graphics_color('bluish green'));
  patch('Faces',R(Utot>1,:),'Vertices',V, ...
  'FaceColor',0.6*[0,158,115]./255);
  patch('Faces',R(Utot>0&Q(2,:,i)>kappa_prol,:),'Vertices',V, ...
  'FaceColor',graphics_color('vermillion'));
  patch('Faces',R(Utot>1&Q(2,:,i)>kappa_prol,:),'Vertices',V, ...
  'FaceColor',0.6*[213,94,0]./255);
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V,...   % necrotic cells
  'FaceColor','black');
  axis(0.4*[-1, 1, -1, 1])
  axis off
  drawnow
end

% -- uncomment to save fig as pdf:
%exportgraphics(f,"tumorfigNN.pdf" , 'resolution', 600)
