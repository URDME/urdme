%% Downsclaled AvascularTumour problem 11x11 with inner 3x3 geometry

% Number of voxels and step size h
Nvoxels = 43;
h = 2/(Nvoxels-1);

% Choose to add idof3
do_idof3 = true;

% Set up mesh and matrixes
[P,E,T,gradquotient] = basic_mesh2(1,Nvoxels);
[L,dM,N,M] = dt_operators(P,T);
Mgamma = assemble_Mgamma(P,T);
[L_orig,M_orig] = assema(P,T,1,1,0);
[V,R] = mesh2dual(P,E,T,'voronoi');
Mgamma_dM = Mgamma./dM;
N_Mgamma = (Mgamma ~= 0) - speye(size(Mgamma));
neigh = full(sum(N,2));

% Initial population that creates 3x3 inner geometry
IC = 5; % Choose initial condition (1,2,3,4,5,6)
R1 = 0.3; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)
Pr = setInitialCondition(IC,R1,R2,P,Nvoxels);
% Pr = circshift(Pr,-1);
VPr = (Pr ~= 0);
% Pr(Pr == -1) = 0; % Set all dead cells to 0 to create internal boundary cells (idof2)

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

% Get all neihgbours to idof1:s that should be added another -1/h^2
% find all neighbours to idof1
[r_ind,c] = find(N(idof1,:));
% only keep those that are un-active (i.e. not Adof)
keep = ismember(c,[adof;idof]); 
r_ind = reshape(r_ind(~keep),[],1); c = reshape(c(~keep),[],1);
% Get actual rows from idof1
r = idof1(r_ind);
% Shift the column values to the opposite position around the diagonal
c_shift = -c + r*2;
% Get new matrix with 1 for all neighbours that should be scaled
actvN_scale = fsparse(r, c_shift, 1, size(N));
% Multiply all rows with three unactive neighbours (where the active
% neighbour should be scaled by 3)
ind_3neighs = (sum(actvN_scale,2) == 3);
actvN_scale(ind_3neighs,:) = actvN_scale(ind_3neighs,:)*3;

% Determine also a local enumeration, eg. [1 2 3
% ... numel(Adof)].
Adof = [adof;idof];
 Adof_ = (1:numel(Adof))';
[adof_,sdof_,idof_,idof1_,idof2_,bdof_m_] = ...
      map(Adof_,Adof,adof,sdof,idof,idof1,idof2,bdof_m);

% Plot how the start looks like
plotInitialPic;

% Step through different values of alpha
alpha = 10.^(-4:0.2:4);
alpha_inv_vec = [1./alpha,0];
ymax = 0;
zmax = 0;
for alpha_inv = alpha_inv_vec

    %%% LHS
    % Get local L for all active dofs
    LaX = L(Adof,Adof);
    % Create matrix with ones on the idof1 diagonal
    Lai = fsparse(idof1_,idof1_,1,size(LaX));
    a_Lai = speye(size(Lai)) - Lai;
    
    LaX = LaX - Lai.*LaX;
    LaX = LaX - diag(sum(Lai*LaX,2));
    
    % Get local Mgamma for all active dofs
    Mgamma_b = Mgamma_dM(Adof,Adof);

    % Get only idof1 part of Mgamma_b and then the diagonal
    Mgamma_b_toAdd = Lai*Mgamma_b*Lai - Mgamma_b.*Lai;
    Mgamma_b_toAdd = Mgamma_b_toAdd + diag(2*sum(Mgamma_b_toAdd,2));
    
    % Adds a -1/h^2 to the L elements that should have one
%     to_add = actvN_scale(Adof,Adof);
%     LaX = LaX + (LaX.*to_add);
    lhs = LaX - a_Lai*LaX*Lai;

    %%% RHS   
    rhs = fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]);       
%     rhs = fsparse(Adof_,1,(cos(2*pi*P(1,Adof)).*sin(2*pi*P(2,Adof)))',[size(LaX,1) 1]);
%     rhs = rhs - mean(rhs); 

    %%% SOLVE
    X_Pr = full(lhs \ rhs);
%     assert(sum(X_Pr < 0) == 0, "X_Pr solution gives negative pressures!")

    %%% PLOT PRESSURE DIFF
    fig1 = figure(1);
    fig1.Position = [0 300 500 400];
    [iii,jjj_] = find(N(idof1,adof)); % neighbours...
    grad = max(X_Pr(jjj_) - X_Pr(idof1_(iii)),0);
    bar(grad);
    % Fix data tip labels to show pressure between which nodes
    labels = compose('%d -> %d', [idof1_(iii),jjj_]);
    dcm_obj = datacursormode(fig1);
    set(dcm_obj,'UpdateFcn',{@datacursor,labels})
    % set ylim and title
    xlabel('Voxel connections');
    ylabel('Pr difference');
    if ymax < max(grad)
        ymax = max(grad);
    end
    ylim([0 ymax]);
    grid on;
    title(sprintf('sum bars = %d \n', sum(grad)));
    
    %%% PLOT SURF
    fig2 = figure(2);
    fig2.Position = [500 300 500 400];
    ax = gca;
    Pr_ = zeros(size(Pr)); Pr_(Adof) = X_Pr;
    Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
    Pr_reshape(Pr_reshape == 0) = NaN;
    [x_Pr,y_Pr] = meshgrid(linspace(-1,1,size(Pr_reshape,1)));
    surf(ax,x_Pr,y_Pr,Pr_reshape,'FaceAlpha', 1.0);
    if zmax < max(X_Pr)
        zmax = max(X_Pr);
    end
%     zlim([0 zmax])
    hold off;
    grid on;
%     colormap('jet');
    condi = rcond(full((lhs)));
    title(sprintf('alpha = %e \n rcond = %e', 1/alpha_inv,condi));
    view(3);
    
    % Print condition number of LHS with rcond
    fprintf('alpha = %e \t rcond(lhs) = %e  \n',1/alpha_inv,condi);
    if condi < 1e-15
        pause
        disp('Not very good this LHS huh?')
    end
    pause(0.1)
end

%% Some display of h's value
disp('----------------------------------------');
fprintf('h = %s \t h^2 = %s \t 1/h = %s \t 1/h^2 = %s \n', ...
    strtrim(rats(h)),strtrim(rats(h^2)), strtrim(rats(1/h)),strtrim(rats(1/h^2)));
disp('----------------------------------------');