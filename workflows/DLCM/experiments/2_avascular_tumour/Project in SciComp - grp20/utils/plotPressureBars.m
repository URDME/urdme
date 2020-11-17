% Function that plots the pressure over Adof with a 3D bar plot
function plotPressureBars(P,U,Pr,adof,adof_,idof,idof_,Nvoxels, downSamp)
    
    % Get the pressure for the active dofs
    Pr_ = full(U); Pr_(adof) = Pr(adof_);
    
    % Reshape the pressure matrix to match the mesh environment
    Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
    
    % If specified, downsample grid to minimize plotting time and comp.
    % efffort
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
    map_start = [0,0,1];
    map_stop = [0.7,0.7,1];
    xx = linspace(0,1,10);
    map_matrix = map_start' + xx.*(map_stop' - map_start');
    mymap = map_matrix';
    caxis([1.1*min(min(Pr_reshape)) 0]);
    colormap(mymap)
    
    title('Pressure in adof(green/orange) and idof(blue)')
    hold off;
    view(3)
end