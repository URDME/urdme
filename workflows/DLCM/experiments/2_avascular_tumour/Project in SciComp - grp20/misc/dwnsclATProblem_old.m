%% Downsclaled AvascularTumour problem 11x11 with inner 3x3 geometry

% Number of voxels and step size h
Nvoxels = 61;
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
IC = 2; % Choose initial condition (1,2,3,4,5,6)
% R1 = 2*h*0.99; % Radius of whole initial tumour
% R2 = h*0.99; % Radius of inner initial setup (doubly occupied, dead etc.)
R1 = 0.5; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)
Pr = setInitialCondition(IC,R1,R2,P,Nvoxels);
% P_inds = find((P(1,:) > -1 & P(1,:) < 1).*(P(2,:) > -1 & P(2,:) < 1));
% Pr = fsparse(P_inds(:),1, ones(size(P_inds(:))), [Nvoxels^2 1]); 
% 
% Pr(1800) = 2;
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


% Determine also a local enumeration, eg. [1 2 3
% ... numel(Adof)].
Adof = [adof;idof];
 Adof_ = (1:numel(Adof))';  
[adof_,sdof_,idof_,idof1_,idof2_,bdof_m_] = ...
      map(Adof_,Adof,adof,sdof,idof,idof1,idof2,bdof_m);


% Plot how the start looks like
shitForMakingInitialPic;

% Get all neihgbours to idof1:s that should be added another -1/h^2
LAI = fsparse(idof1,idof1,1,size(N));
[r,c] = find(LAI*N);
keep = ismember(c,Adof);
r = reshape(r(keep),[],1); c = reshape(c(keep),[],1);
neighs = fsparse(r,c,1,size(N));
unactive_neighs = LAI*N-neighs;
[r_un, c_un] = find(unactive_neighs);
c_un_shift = -c_un + r_un*2;
unactive_neighs_shift = fsparse(r_un, c_un_shift, 1, size(N));
to_add = (unactive_neighs_shift-unactive_neighs)>0;
to_add = to_add(Adof,Adof);

% Step through different values of alpha
alpha = [1e+1, 4/h^2, 1e+5, 1e+7];
alpha_inv = 1./alpha;
i = 1;
ymax = 0;
zmax = 0;
g = 0;
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
    Mgamma_b_toAdd = Lai*Mgamma_b*Lai - Mgamma_b.*Lai;

%     Mgamma_b_toAdd = Mgamma_b + diag(2*sum(Mgamma_b,2)); % PROBABLY THIS
    Mgamma_b_toAdd = diag(2*sum(Mgamma_b_toAdd,2));  % BUT MAYBE THIS

    % Count the number neighs to all Adofs
    neighs_LaX = sum(LaX~=0,2)-1;
    % Create a scaling matrix which scales the hat functions depending on
    % active dofs
    scale_LaX = 4*(ones(size(neighs_LaX)) - neighs_LaX/4);
    scale_LaX_1quart = scale_LaX == 1;
    scale_LaX_halfs = scale_LaX == 2;
    scale_LaX_3quart = scale_LaX == 3;
    
%     Mgamma_b_toAdd = diag(scale_LaX)*2*h - diag(scale_LaX).*scale_LaX_3quart*2/3*h;
    
    % Adds a -1/h^2 to the L elements that should have one NEEDED!
    LaX = LaX + (LaX.*to_add) + 2*LaX.*(to_add.*scale_LaX_3quart);% - a_Lai*LaX*Lai;
    % Scale the boundary rows such that no non-diag L element is less than -1/h^2 MAYBE NOT NEEDED!
%     lhs = LaX - 3/4*Lai*LaX.*scale_LaX_1quart - 3/4*Lai*LaX.*scale_LaX_halfs - 3/4*Lai*LaX.*scale_LaX_3quart + a_inv*Mgamma_b_toAdd;
%     lhs = LaX - a_Lai*LaX*(Lai.*scale_LaX)*1/4 + a_inv*Mgamma_b_toAdd;
    lhs = LaX + a_inv.*Mgamma_b_toAdd;
    % Create final lhs
%     lhs = LaX - 3/4*a_Lai*LaX - 3/4*Lai*LaX.*scale_LaX_1quart - 3/4*Lai*LaX.*scale_LaX_halfs - 3/4*Lai*LaX.*scale_LaX_3quart + a_inv*Mgamma_b_toAdd;

    % Put Dirichlet BC on inner boundary idof2
%     Lai2 = fsparse(idof2_,idof2_,1,size(lhs));
%     lhs = lhs - Lai2*lhs + Lai2;

    %%% RHS   
%     rhs = fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]) + ...
%         fsparse(idof1_,1,ones(size(idof1_))*g*a_inv/h,[size(LaX,1) 1]);

%     rhs = fsparse(Adof_,1,(cos(2*pi*P(1,Adof)).*sin(2*pi*P(2,Adof)))',[size(LaX,1) 1]);
    rhs = fsparse(sdof_,1,1./dM(sdof),[size(LaX,1) 1]);

    %%% SOLVE
    X_Pr = full((lhs'*lhs) \ (lhs'*rhs));
%     X_Pr = normalize(X_Pr,'range');
    assert(sum(X_Pr < 0) == 0 || sum(X_Pr(X_Pr < 0)) < 1e-4, "X_Pr solution gives negative pressures!")

    %%% PLOT
    figure(1);
    subplot(2,2,i);   
    
    [iii,jjj_] = find(N(idof1,adof)); % neighbours...
    grad = max(X_Pr(jjj_) - X_Pr(idof1_(iii)),0);
    bar(grad);
%     xticks(1:length(iii))
% %     xticklabels(compose('%d\\newline->\\newline%d', [jjj_,idof1_(iii)]));
%     xt = get(gca, 'XTick');
%     labels = compose('%d\\newline->\\newline%d', [jjj_,idof1_(iii)]);
%     text(xt, grad, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',9)
%     set(gca,'xticklabel',{[]})
    ylabel('Pr diff');
    if ymax < max(grad)
        ymax = max(grad);
    end
    ylim([0 ymax]);
    grid on;
    title(sprintf('sum bars = %d \n', sum(grad)));
    
    
    %%% PLOT SURF
    fig2 = figure(2);
    subplot(2,2,i);
    ax = gca;
%     plotPressureBars2(fig2,Nvoxels,h,Pr,X_Pr,adof,adof_,idof1,idof1_)
    Pr_ = zeros(size(Pr)); Pr_(Adof) = X_Pr;
    Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
    Pr_reshape(Pr_reshape == 0) = NaN;
    % Generate grid points  
    [x_Pr,y_Pr] = meshgrid(linspace(-1,1,size(Pr_reshape,1)));
    surf(ax,x_Pr,y_Pr,Pr_reshape,'FaceAlpha', 1.0);
%     hold on;
%     Pr_ = zeros(size(Pr)); Pr_(idof1(iii)) = grad;
%     Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
    % Generate grid points  
%     scatterbar3(x_Pr,y_Pr,Pr_reshape,h/2);
    if zmax < max(X_Pr)
        zmax = max(X_Pr);
    end
%     zlim([0 zmax])
%     axis([-1 1 -1 1 0 0.2]);
    hold off;
    grid on;
    title(sprintf('alpha = %d', 1/a_inv));
%     colormap('jet')
    view(3);
    
%     %%% PLOT GRADIENT
%     figure(3);
%     subplot(2,2,i);
%     Pr_ = full(Pr); Pr_(Adof) = X_Pr;
%     z = reshape(Pr_, Nvoxels, Nvoxels);
%     [x,y] = meshgrid(linspace(-1,1,size(z,1)));
%     [px,py] = gradient(z);
%     contour(x,y,z)
%     hold on
%     quiver(x,y,px,py)
%     hold off
%     grid on;

    % Step figure index
    i = i + 1;
end
% Create figure handles and set names
fig1 = figure(1);
sgtitle(fig1, 'Diff in Pr in boundary')
fig1.Name = 'diffPr-plot';
fig2 = figure(2);
sgtitle(fig2, 'Pr in adof and idof')
fig2.Name = 'surf-plot';
% fig3 = figure(3);
% sgtitle(fig3, 'Gradient for Pr')
% fig3.Name = 'gradient-plot';

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