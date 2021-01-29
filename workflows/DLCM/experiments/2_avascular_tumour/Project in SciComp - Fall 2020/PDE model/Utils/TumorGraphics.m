% Can be run after AvascularTumour for additional visuals 

snapshotfigure = 0; %creates a gif and 5 snapshot figures
doffigure = 1;      %visualizes the sdof's and bdof's as well as the concentration in other voxels
deadfigure = 1;     %visualizes the dead cells in every voxel 
oxyfigure = 1;      %visualizes the concentration oxygen in the grid

%% Snapshot figure
if snapshotfigure == 1
    % create a GIF animation
    Tumour = struct('cdata',{},'colormap',{});
    fig=figure(1);
    clf,

    Umat=full(cell2mat(Usave));
    colorbar;
    caxis([0 max(max(Umat))])
    colorlabel('Concentration of cells, U')

    snapshot = 0;    %Save 5 snapshots 

    for i = 1:numel(Usave)        
        % background
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');

        hold on,
        axis([-1 1 -1 1]); axis square, axis off

        % colour living voxels after concentration level
        ii = find(Usave{i}>0);
        c = Umat(ii,i);
        patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c, ... 
            'FaceColor','flat');     

        % colour (fully) dead voxels black
        ii = find(Usave{i} == 0 & Udsave{i} > 0);
        p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
        %legend(p_dead,'dead')

        title(sprintf('Time = %d',tspan(i)));
        drawnow;
        Tumour(i) = getframe(gcf);

        %Save 5 snapshots of the tumor progression
        if i~= [1 ceil([0.24 0.49 0.74 1]*numel(Usave))]
        elseif i==1
            filename = 'T=1.png';
            print(fig,filename,'-painters','-dpng');
        else
            ii=ceil(i*(Tend/numel(Usave)));
            filename = ['T=' num2str(ii) '.png'];
            print(fig,filename,'-painters','-dpng');
        end


        % saves the GIF
        movie2gif(Tumour,{Tumour([1:2 end]).cdata},'Tumour.gif', ...
                 'delaytime',0.1,'loopcount',0);

    end
end

%% dof figure
if doffigure==1
    % create a GIF animation
    Mdof = struct('cdata',{},'colormap',{});
    figure(2), clf,

    Umat=full(cell2mat(Usave));
    colorbar
    caxis([0 max(max(Umat))])
    colorlabel('Concentration of cells, u')
    for i = 1:numel(Usave)

        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off

        ii = find(Usave{i}>0);
        c = Umat(ii,i);
        patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
        
        % colour the different dof types
        p_bdof = patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
            'FaceColor','cyan'); 

        p_sdof = patch('Faces',R(sdofsave{i},:),'Vertices',V, ...
            'FaceColor','magenta'); 

        % colour (fully) dead voxels black
        ii = find(Usave{i} == 0 & Udsave{i} >0);
        p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
        legend([p_bdof,p_sdof,p_dead],'bdof','sdof','dead')

        title(sprintf('Time = %d',tspan(i)));

        drawnow;
        Mdof(i) = getframe(gcf);
    end

    % saves the GIF
    movie2gif(Mdof,{Mdof([1:2 end]).cdata},'TumourMdof.gif', ...
              'delaytime',0.1,'loopcount',0);
end
    
    
%% oxygen figure
if oxyfigure==1
    % create a GIF animation
    Mnormal = struct('cdata',{},'colormap',{});
    figure(3), clf,

    Oxymat=full(cell2mat(Oxysave));
    colorbar
    caxis([min(min(Oxymat)) max(max(Oxymat))])
    colorlabel('Concentration of oxygen, c')
    for i = 1:numel(Usave)

        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off

        % Colour of oxygen concentration
        ii = find(Oxysave{i}>0);
        c = Oxymat(ii,i);
        patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');     

        title(sprintf('Time = %d',tspan(i)));
        drawnow;
        Mnormal(i) = getframe(gcf);
    end
    % saves the GIF
    movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'Oxynormal.gif', ...
              'delaytime',0.1,'loopcount',0);
end

%% Dead figure
if deadfigure==1
    Mdead = struct('cdata',{},'colormap',{});
    figure(4), clf,

    Udmat=cell2mat(Udsave);
    cmat = full(Udmat/max(max(Udmat)));
    caxis([1 2])
    colorbar;
    colormap 'gray'

    for i = 1:numel(Udsave)
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off

        % colour (fully) dead voxels black
        ii = find(Udsave{i}>0);
        c = cmat(ii,i);
        patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    

        title(sprintf('Time = %d',tspan(i)));
        drawnow;
        Mdead(i) = getframe(gcf);
    end
    % saves the GIF
    movie2gif(Mdead,{Mdead([1:2 end]).cdata},'TumourMDead.gif', ...
              'delaytime',0.1,'loopcount',0);
end

