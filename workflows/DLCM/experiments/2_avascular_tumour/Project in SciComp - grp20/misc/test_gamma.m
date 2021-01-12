% ---------------------------------------------------------------------- %
% Test down-scaled problem
% ---------------------------------------------------------------------- %

Nvoxels = 11;
h = 2/(Nvoxels-1);

[p,e,t,gradquotient] = basic_mesh(1,Nvoxels);
[L,dM,N,M] = dt_operators(p,t);
Mgamma = assemble_Mgamma2(p,t);
Mgamma = Mgamma./dM;
neigh = full(sum(N,2));
N_Mgamma = (Mgamma ~= 0) - speye(size(Mgamma));
[V,R] = mesh2dual(p,e,t,'voronoi');

% Initial population that creates 3x3 inner geometry
IC = 2; % Choose initial condition (1,2,3,4,5,6)
R1 = 0.3; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)
Pr = setInitialCondition(IC,R1,R2,p,Nvoxels);
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

do_idof3 = 1;
if do_idof3
    idof1 = sort([idof1;idof3]);
    idof = sort([idof;idof3]);
end

% Determine also a local enumeration, eg. [1 2 3]
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
alpha = [1e-2, 1e-1, 1e0, 1e2];
alpha_inv = 1./alpha;
i = 1;
for a_inv = alpha_inv

    %%% LHS
    % Get local L for all active dofs
    LaX = L(Adof,Adof);
    
    % Get current Mgamma for all active dofs
    Mgamma_b = Mgamma(Adof,Adof);
    
    % Create matrix with ones on the idof1 diagonal
    Lai = fsparse(idof1_,idof1_,1,size(LaX));

    % Get only idof1 part of Mgamma_b
    Mgamma_b = Lai*Mgamma_b*Lai;
    
    % Add edge mass contribution to idof nodes
    Mgamma_toAdd = diag(sum(Mgamma_b,2));

%     % Find adofs which have neighbouring idofs
%     b_adof = LaX(idof_,adof_);
%     
%     % Find adof index and occurance of idof neighbours
%     [~,j] = find(b_adof);
%     [N_idof_neigh, adof_ind] = groupcounts(j);
%     
%     % Create vector that determines how many idof neighbours one adof has
%     adof_scale = zeros(length(adof_),1);
%     adof_scale(adof_ind) = N_idof_neigh;
%     
%     % Create a scaling matrix which scales the hat functions depending on
%     % the number of idof neighbours
%     full_neigh = fsparse(adof_,adof_, 4, size(LaX));
%     scale_neigh = fsparse(adof_,adof_, adof_scale, size(LaX));
%     
%     % Scale the Laplacian
%     LaX = (1/4).*(full_neigh - scale_neigh)*LaX;

    % Find how many connections idof1 have
    [~,j] = find(LaX(numel(adof_)+1:end, numel(adof_)+1:end));
    j_ = j + numel(adof);
    
    % Find adof index and occurance of idof neighbours
    [N_neigh, idof_ind] = groupcounts(j_);
    
    % Create vector that determines how many neighbours one idof has
    idof_scale = zeros(length(idof_),1);
    idof_scale(idof_ind-numel(adof_)) = (1-(1/4)*(4-N_neigh));
    
    % Scale Laplacian
    LaX(idof_,idof_) = repmat(idof_scale,1,numel(idof_))'.*LaX(idof_, idof_);
    
    % idof neighbours to adof part
    LaX(numel(adof_)+1:end, 1:numel(adof_)) = repmat(idof_scale,1,numel(adof_)).*LaX(numel(adof_)+1:end, 1:numel(adof_)); 
%         LaX(numel(adof_)+1:end, 1:numel(adof_)) = 0*LaX(numel(adof_)+1:end, 1:numel(adof_));
        
%     % adof neighbours to idof part
%     LaX(1:numel(adof_), numel(adof_)+1:end) = repmat(idof_scale,1,numel(adof_))'.*LaX(1:numel(adof_), numel(adof_)+1:end);
%       
     LaX(1:numel(adof_), numel(adof_)+1:end) = 0.*LaX(1:numel(adof_), numel(adof_)+1:end);

    %%% LHS
    lhs = LaX + a_inv*Mgamma_toAdd;
    
    %%% RHS   
    rhs = ones(size(LaX,1),1) + fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]);

    %%% SOLVE
    X_Pr = full(lhs \ rhs);

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



% % ---------------------------------------------------------------------- %
% % Test matrix down-scaled problem
% % ---------------------------------------------------------------------- %
% 
% % Parameters 
% Nvoxels = 2;
% alpha = 10.^(-3:5);
% alpha_inv = 1./alpha;
% 
% % Create mesh and assemble matrices
% [P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
% [V,R] = mesh2dual(P,E,T,'voronoi');
% 
% % Multiply L by inverse of M
% [L0,M] = assema(P,T,1,1,0);
% 
% % Invert L using lumped mass matrix
% % the (lumped) mass matrix dM gives the element volume
% dM = full(sum(M,2));
% ndofs = size(dM,1);
% 
% % explicitly invert the lumped mass matrix and filter the diffusion matrix
% [i,j,s] = find(L0);
% s = s./dM(i);
% keep = find(i ~= j); % (removes only the diagonal)
% i = i(keep); j = j(keep); s = s(keep);
% 
% % rebuild L, ensuring that the diagonal equals minus the sum of the
% % off-diagonal elements
% LL = sparse(i,j,s,ndofs,ndofs);
% LL = LL+sparse(1:ndofs,1:ndofs,-full(sum(LL,2)));
% 
% N = sparse(i,j,1,ndofs,ndofs);
% 
% % Assemble boundary matrix and filter it as the
% Mgamma = assemble_Mgamma(P,T);
% Mgamma_dM = Mgamma./dM;
% 
% % Create DOFs
% adof = 1:floor(Nvoxels^2/2); % Active voxels (source voxels)
% idof = adof(end)+1:size(LL,1); 
% 
% % Laplacian adjustment
% Lai = fsparse(idof,idof,1,size(LL));
% 
% % Investigate matrix multiplication for different values of alpha
% RHS = ones(size(LL,1),1) + full(fsparse(adof',1,1./dM(adof),[size(LL,1) 1]));
% for a_inv = 1e6%alpha_inv
%     
%     % Case one: All matrices being inverted at the same time
%     LHS1 = (M \ (L0 - Lai*L0 + a_inv*Lai*Mgamma));
%     X1 = LHS1 \ RHS;
%     
%     % Case two: L and M inverted in dt_operators.m
%     LHS2 = LL - Lai*LL + a_inv*Lai*Mgamma_dM;
%     X2 = LHS2 \ RHS;
%     
% end

% ---------------------------------------------------------------------- %
% Test M_Gamma function
% ---------------------------------------------------------------------- %

% [p,e,t] = simple_mesh(3);
% 
% % Mgamma-function
% np = size(p,2); % Number of nodes 
% nt = size(t,2); % Number of boundary edges
% Mgamma = zeros(np,np); % Allocate mass matrix 
% 
% for N = 1:np
%     
%     % Find all edges that are connected to N
%     % Find column number (triangle) in which the node N is part of 
%     [~, col] = find(t(1:3,:)==N);
%     
%     % Find other nodes in these triangles
%     tri = t(1:3, col);
%     N_neigh = unique(tri(tri~=N));
% 
%     % Calculate edge lengths from N to its neighbours
%     len = sqrt((p(1,N) - p(1,N_neigh)).^2 + ((p(2,N) - p(2,N_neigh)).^2));
%     
%     Mgamma(N,N_neigh) = (1/6).*len;
%     Mgamma(N,N) = (2/6);
%     
% end
