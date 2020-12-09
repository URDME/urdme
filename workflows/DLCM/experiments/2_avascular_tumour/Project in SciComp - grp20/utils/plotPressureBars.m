%% Plots the pressure over Adof with a 3D bar plot
    
% Get the pressure for the active dofs
Pr_ = full(U); Pr_(adof) = Pr(adof_);

% Reshape the pressure matrix to match the mesh environment
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);

% If specified, downsample grid to minimize plotting time and comp.
% efffort
downSamp = 1;
if downSamp < 1
    Pr_reshape = imresize(Pr_reshape, downSamp, 'bilinear');
end

% Generate grid points  
[x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,size(Pr_reshape,1)));

% Set width of bars
width = sqrt((P(1,1) - P(1,2))^2 + (P(2,1) - P(2,2))^2);

% Patch and plot bars
hold on;
scatterbar3(x_Pr_,y_Pr_,Pr_reshape,width);
grid on;

% Set colormap for the active dofs
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
% caxis([1.1*min(Pr_(ii)) 1.1*max(Pr_(ii))]);
caxis([mean(Pr_(ii))-2*std(Pr_(ii)) mean(Pr_(ii))+2*std(Pr_(ii))]);
% caxis([0 3]);
colormap(mymap)

% Freeze this colormap (in order to apply another one to the boundary
% dof)
freezeColors;

% Get the pressure for the boundary dofs
Pr_(adof) = 0;
Pr_(idof) = Pr(idof_);

% Reshape the pressure matrix to match the mesh environment
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);

% If specified, downsample grid to minimize plotting time and comp.
% efffort
if downSamp < 1
    Pr_reshape = imresize(Pr_reshape, downSamp, 'bilinear');
end

% Patch and plot bars
scatterbar3(x_Pr_,y_Pr_,Pr_reshape,width);

% Set colormap for the boundary dofs
% caxis([1.1*min(min(Pr_(ii))) 1.1*max(max(Pr_(ii)))]);
caxis([mean(Pr_(ii))-2*std(Pr_(ii)) mean(Pr_(ii))+2*std(Pr_(ii))]);
% caxis([-1 0.5]);
colormap('cool')

title(sprintf('Time = %d, Ncells = %d \n alpha = %d' , ...
            tspan(end),full(sum(abs(Usave{end}))), alpha));
hold off;
view(3)