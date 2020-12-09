%% Downsclaled AvascularTumour problem 11x11 with inner 3x3 geometry

% Number of voxels and step size h
Nvoxels = 11;
h = 2/(Nvoxels-1);

% Set up mesh and matrixes
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[L,dM,N,M] = dt_operators(P,T);
[L_orig,M_orig] = assema(P,T,1,1,0);
[V,R] = mesh2dual(P,E,T,'voronoi');
Mgamma = assemble_Mgamma(P,T);
Mgamma_dM = Mgamma./dM;
neigh = full(sum(N,2));

% Initial population that creates 3x3 inner geometry
IC = 2; % Choose initial condition (1,2,3,4,5,6)
R1 = 0.35; % Radius of whole initial tumour
R2 = 0.10; % Radius of inner initial setup (doubly occupied, dead etc.)
Pr = setInitialCondition(IC,R1,R2,P,Nvoxels);
VPr = (Pr ~= 0);
% Pr(Pr == -1) = 0; % Set all dead cells to 0 to create internal boundary cells (idof2)

adof = find(Pr); % all filled voxels
bdof_m = find(N*(Pr ~= 0) < neigh & abs(Pr) >= 1); % single occupied cells next to external boundary
sdof = find(Pr > 1); % voxels with 2 cells
Idof = (N*(Pr ~= 0) > 0 & Pr == 0); % empty voxels touching occupied ones
idof1 = find(Idof & ~VPr); % "external" BC1
idof2 = find(Idof & VPr); % "internal" BC2
idof = find(Idof);

% Determine also a local enumeration, eg. [1 2 3
% ... numel(Adof)].
Adof = [adof;idof];
 Adof_ = (1:numel(Adof))';  
[adof_,sdof_,idof_,idof1_,idof2_,bdof_m_] = ...
      map(Adof_,Adof,adof,sdof,idof,idof1,idof2,bdof_m);

% Plot how the start looks like
figure(12);
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
ii = find(Pr == 1);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'));
ii = find(Pr == 2);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('vermillion'));
ii = find(Pr == -1);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0]);
ii = idof1;
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0,0,1]);
title('Start out structure')
drawnow;

% Step through different values of alpha
alpha = [1e-2, 1e-1, 1e+1, 1e+2];
alpha_inv = 1./alpha;
i = 1;
for a_inv = alpha_inv

    %%% LHS
    % Get local L for all active dofs
    LaX = L(Adof,Adof);
    % Create matrix with ones on the idof1 diagonal
    Lai = fsparse(idof1_,idof1_,1,size(LaX));
    % Get local Mgamma for all active dofs
    Mgamma_b = Mgamma_dM(Adof,Adof);
    % Get only idof1 part of Mgamma_b
    Mgamma_b_toAdd = Lai*Mgamma_b*Lai;
    % Count the number neighs to all Adofs
    neighs_LaX = sum(LaX~=0,2)-1;
    % Create a scaling matrix which scales the hat functions depending on
    % active dofs
    scale_LaX = fsparse(diag(ones(size(neighs_LaX)) - neighs_LaX./4,0));
    
    % Get lhs
    lhs = LaX - Lai*LaX*scale_LaX + a_inv*Mgamma_b_toAdd;
    
    % Put Dirichlet BC on inner boundary idof2
    Lai2 = fsparse(idof2_,idof2_,1,size(lhs));
    lhs = lhs - Lai2*lhs + Lai2;

    %%% RHS   
    rhs = fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]);
    
    %%% SOLVE
    X_1 = full(lhs \ rhs);
%     X_1 = normalize(X_1,'range');

    %%% PLOT
    figure(1);
    subplot(2,2,i);   
    plot(bdof_m, X_1(bdof_m_),'.-', 'Displayname', 'bdof_m');
    hold on;
    plot(idof, X_1(idof_),'.-', 'Displayname', 'idof1');
    title(sprintf('alpha = %d', 1/a_inv));
    grid on;

    %%% PLOT SURF
    fig2 = figure(2);
    subplot(2,2,i);
    plotPressureBars2(fig2,Nvoxels,h,Pr,X_1,adof,adof_,idof,idof_)
    axis([-1 1 -1 1]);
    grid on;
    title(sprintf('alpha = %d', 1/a_inv));
    view(3);

    % Step figure index
    i = i + 1;
end

% Plot scaling of matrixes used in problem
row = idof1(1);
format rational
disp('----------------------------------------');
fprintf('NVOXELS = %d --> size(L) = %dx%d\n', Nvoxels,Nvoxels^2,Nvoxels^2);
fprintf('ROW IN MATRIXES: %d\n', row);
fprintf('h = %s \t h^2 = %s \t 1/h = %s \t sqrt(2)*h = %s \t 2*(2+sqrt(2))*h/6 = %s \n', ...
    strtrim(rats(h)),strtrim(rats(h^2)), strtrim(rats(1/h)), ...
    strtrim(rats(sqrt(2*h^2))), strtrim(rats(2*(2+sqrt(2))*h/6)));
fprintf('L \t\t\t= %s\n', strtrim(rats(full(L(row,L(row,:) ~= 0)))));
fprintf('L_orig \t\t= %s\n', strtrim(rats(full(L_orig(row,L_orig(row,:) ~= 0)))));
fprintf('Mgamma \t\t= %s\n', strtrim(rats(full(Mgamma(row,Mgamma(row,:) ~= 0)))));
fprintf('Mgamma_dM \t= %s\n', strtrim(rats(full(Mgamma_dM(row,Mgamma_dM(row,:) ~= 0)))));
fprintf('M \t\t\t= %s\n', strtrim(rats(full(M(row,M(row,:) ~= 0)))));
fprintf('dM \t\t\t= %s\n', strtrim(rats(dM(row))));
disp('----------------------------------------');
format short