%% Downsclaled AvascularTumour problem 11x11 with inner 3x3 geometry

% Number of voxels and step size h
Nvoxels = 11;
h = 2/(Nvoxels-1);

% Choose to add idof3
do_idof3 = true;

% Set up mesh and matrixes
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
% [~,~,T_2,~] = flipped_mesh(Nvoxels);
% T = [T, T_2];
[L,dM,N,M] = dt_operators(P,T);
Mgamma = assemble_Mgamma(P,T);
% Mgamma(Mgamma < 0.0333 & Mgamma > 0) = 197/4179;
% dM = dM/2; M = M/2; Mgamma = Mgamma/2;
[L_orig,M_orig] = assema(P,T,1,1,0);
[V,R] = mesh2dual(P,E,T,'voronoi');
Mgamma_dM = Mgamma./dM;
N_Mgamma = (Mgamma ~= 0) - speye(size(Mgamma));
neigh = full(sum(N,2));

% Initial population that creates 3x3 inner geometry
IC = 2; % Choose initial condition (1,2,3,4,5,6)
% R1 = 2*h*0.99; % Radius of whole initial tumour
% R2 = h*0.99; % Radius of inner initial setup (doubly occupied, dead etc.)
R1 = 0.3; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)
Pr = setInitialCondition(IC,R1,R2,P,Nvoxels);
VPr = (Pr ~= 0);
Pr(Pr == -1) = 0; % Set all dead cells to 0 to create internal boundary cells (idof2)

adof = find(Pr); % all filled voxels
bdof_m = find(N*(Pr ~= 0) < neigh & abs(Pr) == 1); % single occupied cells next to external boundary
sdof = find(Pr > 1); % voxels with 2 cells
Idof = (N*(Pr ~= 0) > 0 & Pr == 0); % empty voxels touching occupied ones
idof1 = find(Idof & ~VPr); % "external" BC1
idof2 = find(Idof & VPr); % "internal" BC2
idof3 = find(~VPr & N_Mgamma*VPr > 0); % boundary around "visited voxels"
idof3 = setdiff(idof3,idof1);
idof = find(Idof);

if do_idof3
    idof1 = sort([idof1;idof3]);
    idof = sort([idof;idof3]);
end


% Determine also a local enumeration, eg. [1 2 3
% ... numel(Adof)].
Adof = [adof;idof];
 Adof_ = (1:numel(Adof))';  
[adof_,sdof_,idof_,idof1_,idof2_,bdof_m_] = ...
      map(Adof_,Adof,adof,sdof,idof,idof1,idof2,bdof_m);


% Plot how the start looks like
figure(99);
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
hold off;
drawnow;

% Step through different values of alpha
alpha = [1e-2, 1e-1, 1e-0, 1e+1];
alpha_inv = 1./alpha;
i = 1;
for a_inv = alpha_inv

    %%% LHS
    % Get local L for all active dofs
    LaX = L(Adof,Adof);
    % Create matrix with ones on the idof1 diagonal
    Lai = fsparse(idof1_,idof1_,1,size(LaX));
    a_Lai = speye(size(Lai)) - Lai;
    
    % Get local Mgamma for all active dofs
    Mgamma_b = Mgamma_dM(Adof,Adof);

    % Get only idof1 part of Mgamma_b
    Mgamma_b = Lai*Mgamma_b*Lai - Mgamma_b.*Lai;

%     Mgamma_b_toAdd = Mgamma_b + Lai*diag(2*sum(Mgamma_b,2)); % PROBABLY THIS
    Mgamma_b_toAdd = diag(sum(Mgamma_b,2));  % BUT MAYBE THIS

    % Count the number neighs to all Adofs
    neighs_LaX = sum(LaX~=0,2)-1;
    % Create a scaling matrix which scales the hat functions depending on
    % active dofs
    scale_LaX = ones(size(neighs_LaX)) - neighs_LaX/4;
    scale_LaX_halfs = scale_LaX == 0.5;
    scale_LaX_1quart = scale_LaX == 0.25;
    scale_LaX_3quart = scale_LaX == 0.75;
    % Adds a -1/h^2 to the L elements that should have one NEEDED!
    LaX = LaX + Lai*LaX*a_Lai.*(scale_LaX_1quart + 3*scale_LaX_3quart) + Lai*(LaX - Lai.*LaX).*scale_LaX_halfs;
%     figure;imagesc(LaX);
    % Scale the boundary rows such that no non-diag L element is less than -1/h^2 MAYBE NOT NEEDED!
    lhs = LaX - 1/2*Lai*LaX.*scale_LaX_1quart - 3/4*Lai*LaX.*scale_LaX_halfs - 3/4*Lai*LaX.*scale_LaX_3quart - a_Lai*LaX*Lai + a_inv*Mgamma_b_toAdd;
%     lhs = LaX - a_Lai*LaX*(Lai.*scale_LaX) + a_inv*Mgamma_b_toAdd; 
%     lhs = LaX + a_inv*Mgamma_b_toAdd;
    % Create final lhs
%     lhs = LaX - 1/2*Lai*LaX.*scale_LaX_1quart - 3/4*Lai*LaX.*scale_LaX_halfs - 3/4*Lai*LaX.*scale_LaX_3quart + a_inv*Mgamma_b_toAdd;

    % Put Dirichlet BC on inner boundary idof2
%     Lai2 = fsparse(idof2_,idof2_,1,size(lhs));
%     lhs = lhs - Lai2*lhs + Lai2;

    %%% RHS   
    rhs = ones(size(LaX,1),1) + fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]);

    %%% SOLVE
    X_Pr = full(lhs \ rhs);
%     X_Pr = normalize(X_Pr,'range');

    %%% PLOT
    figure(1);
    subplot(2,2,i);   
    plot(adof, X_Pr(adof_),'x-', 'Displayname', 'adof');
    hold on;
    plot(idof1, X_Pr(idof1_),'.-', 'Displayname', 'idof1');
    title(sprintf('alpha = %d', 1/a_inv));
%     ylim([0 0.2])
%     hold off;
    grid on;

    %%% PLOT SURF
    fig2 = figure(2);
    subplot(2,2,i);
    plotPressureBars2(fig2,Nvoxels,h,Pr,X_Pr,adof,adof_,idof1,idof1_)

%     axis([-1 1 -1 1 0 0.2]);
    grid on;
    title(sprintf('alpha = %d', 1/a_inv));
    view(3);
    
    %%% PLOT GRADIENT
    figure(3);
    subplot(2,2,i);
    Pr_ = full(Pr); Pr_(Adof) = X_Pr;
    z = reshape(Pr_, Nvoxels, Nvoxels);
    [x,y] = meshgrid(linspace(-1,1,size(z,1)));
    [px,py] = gradient(z);
    contour(x,y,z)
    hold on
    quiver(x,y,px,py)
    hold off
    grid on;

    % Step figure index
    i = i + 1;
end
% Create figure handles and set names
fig1 = figure(1);
sgtitle(fig1, 'Pr in adof(x) and idof(.)')
fig1.Name = 'bdof_m-idof-plot';
fig2 = figure(2);
sgtitle(fig2, 'Pr in adof and idof')
fig2.Name = 'surf-plot';
fig3 = figure(3);
sgtitle(fig3, 'Gradient for Pr')
fig3.Name = 'gradient-plot';

%%
% Plot scaling of matrixes used in problem
row = idof1(1);
format rational
disp('----------------------------------------');
fprintf('NVOXELS = %d --> size(L) = %dx%d\n', Nvoxels,Nvoxels^2,Nvoxels^2);
fprintf('ROW IN MATRIXES: %d\n', row);
fprintf('h = %s \t h^2 = %s \t 1/h = %s \t sqrt(2)*h/6 = %s \t 4*(2+2*sqrt(2))*h/6 = %s \n', ...
    strtrim(rats(h)),strtrim(rats(h^2)), strtrim(rats(1/h)), ...
    strtrim(rats(sqrt(2)*h/6)), strtrim(rats(4*(2+2*sqrt(2))*h/6)));
fprintf('L \t\t\t= %s\n', strtrim(rats(full(L(row,L(row,:) ~= 0)))));
fprintf('L_orig \t\t= %s\n', strtrim(rats(full(L_orig(row,L_orig(row,:) ~= 0)))));
fprintf('Mgamma \t\t= %s\n', strtrim(rats(full(Mgamma(row,Mgamma(row,:) ~= 0)))));
fprintf('Mgamma_dM \t= %s\n', strtrim(rats(full(Mgamma_dM(row,Mgamma_dM(row,:) ~= 0)))));
fprintf('M \t\t\t= %s\n', strtrim(rats(full(M(row,M(row,:) ~= 0)))));
fprintf('dM \t\t\t= %s\n', strtrim(rats(dM(row))));
disp('----------------------------------------');
format short