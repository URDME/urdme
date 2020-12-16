function plotPressureBars2(fig,Nvoxels,h,Pr,X_1,adof,adof_,idof1,idof1_)
fig;
Pr_ = full(Pr); Pr_(adof) = X_1(adof_);
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
% Generate grid points  
[x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,size(Pr_reshape,1)));
% Set width of bars
width = h;
% Patch and plot bars
scatterbar3(x_Pr_,y_Pr_,Pr_reshape,width);

% Set colormap for the active dofs
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
% caxis([mean(Pr_)-2*std(Pr_) mean(Pr_)+2*std(Pr_)]);
colormap(mymap)

% Freeze this colormap (in order to apply another one to the boundary
% dof)
freezeColors;

% Get the pressure for the boundary dofs
Pr_(adof) = 0;
Pr_(idof1) = X_1(idof1_);

% Reshape the pressure matrix to match the mesh environment
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);

% Patch and plot bars
scatterbar3(x_Pr_,y_Pr_,Pr_reshape,h);

% Set colormap for the boundary dofs
% caxis([mean(Pr_)-2*std(Pr_) mean(Pr_)+2*std(Pr_)]);
colormap('cool')
end