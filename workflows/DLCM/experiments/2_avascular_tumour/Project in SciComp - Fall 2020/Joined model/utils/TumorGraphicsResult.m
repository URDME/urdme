resultfigure=0; %Visualize 
doffigure = 0;  %visualizes the dead cells in every voxel 
normalfigure = 1;   %visualizes the concentration in every voxel 
deadfigure = 0; %visualizes the sdof's and bdof's as well as the concentration in other voxels
oxyfigure = 0;
%%
if resultfigure==1
% create a GIF animation
Mres = struct('cdata',{},'colormap',{});
figure(25), clf,

Umat=full(cell2mat(Usave));
%colorbar
%caxis([0 max(max(Umat))])
%colorlabel('Concentration of cells, U')
for i = 1:numel(Usave)
    
%     patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%         'EdgeColor','none');
%     hold on,
%     axis([-1 1 -1 1]); axis square, axis off
%     
%     ii = find(Usave{i}>0);
%     c = Umat(ii,i);
%     patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');     
% 
%     ii = find(Usave{i} == 0 & Udsave{i} > 0);
%     p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
%     legend(p_dead,'dead')
%     
    
%%%%%%%%%%%%%%%%%%%%%%%%%%    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = Umat(ii,i);
    
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','g');    
    p_bdof = patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
         'FaceColor','cyan'); 
%     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     p_bdof = patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
%         'FaceColor','cyan'); 
%     
%         
%     
    ii = find(Usave{i}>1);
    c = Umat(ii,i);

    p_sdof = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor','red'); 

    ii = find(Usave{i} == 0 & Udsave{i} >0);
    p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    %legend([p_bdof,p_sdof,p_sdofb,p_dead],'bdof','sdof','sdofb','dead')
    legend([p_bdof,p_sdof,p_dead],'bdof','sdof','dead')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
    Mnormal(i) = getframe(gcf);
    
  %Save 5 snapshots of the tumor progression
  if i~= [1 ceil([0.24 0.49 0.74 1]*numel(Usave))]
  
  elseif i==1
      filename = 'T=1.png';
      saveas(gcf,filename) 
  else
      ii=ceil(i*(Tend/numel(Usave)))
      filename = ['T=' num2str(ii) '.png'];
      saveas(gcf,filename)     
  end
  
end
% saves the GIF
movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'TumourMres.gif', ...
          'delaytime',0.1,'loopcount',0);
end

%% 
if doffigure==1
% create a GIF animation
Mdof = struct('cdata',{},'colormap',{});
figure(1), clf,

Umat=full(cell2mat(Usave));
% cmat = full(Umat/max(max(Umat)));
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
        
%     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
    
    p_bdof = patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
        'FaceColor','cyan'); 
    
        
    %sdof_b
%     p_sdofb = patch('Faces',R(sdofbsave{i},:),'Vertices',V, ...
%         'FaceColor','red'); 
%     
    p_sdof = patch('Faces',R(sdofsave{i},:),'Vertices',V, ...
        'FaceColor','magenta'); 

    ii = find(Usave{i} == 0 & Udsave{i} >0);
    p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    %legend([p_bdof,p_sdof,p_sdofb,p_dead],'bdof','sdof','sdofb','dead')
    legend([p_bdof,p_sdof,p_dead],'bdof','sdof','dead')

    %title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    title(sprintf('Time = %d',tspan(i)));

    drawnow;
    Mdof(i) = getframe(gcf);
end

% saves the GIF
movie2gif(Mdof,{Mdof([1:2 end]).cdata},'TumourMdof.gif', ...
          'delaytime',0.1,'loopcount',0);
end

%%
if oxyfigure==1
% create a GIF animation
Mnormal = struct('cdata',{},'colormap',{});
figure(20), clf,

Oxymat=full(cell2mat(Oxysave));
colorbar
caxis([min(min(Oxymat)) max(max(Oxymat))])
colorlabel('Concentration of oxygen, c')
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Oxysave{i}>0);
    c = Oxymat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');     
% 
%     ii = find(Usave{i} == 0 & Udsave{i} > 0);
%     p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
%     legend(p_dead,'dead')
    
    %title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Oxysave{i})))));
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
    Mnormal(i) = getframe(gcf);
end
% saves the GIF
movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'Oxynormal.gif', ...
          'delaytime',0.1,'loopcount',0);
end
%% BEST PLOT THERE IS!!!
figure(11), clf
set(gca,'color','none');
Umat=full(cell2mat(Usave));
colorbar('southoutside')
caxis([0 1.2])
% caxis([0 max(max(Umat))])
colorlabel('Concentration of cells, U')
axis([-1 1 -1 1]); axis square, axis off
hold on
for i = 1:numel(Usave)
%     patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%         'EdgeColor','none');
    ii = find(Usave{i}>0);
    c = Umat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat','Edgecolor','none');     

    ii = find(Usave{i} == 0 & Udsave{i} > 0);
    p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
end
hold off;
%%
if deadfigure==1
Mdead = struct('cdata',{},'colormap',{});
figure(10), clf,

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
    
    ii = find(Udsave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    

    %title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
    Mdead(i) = getframe(gcf);
end
% saves the GIF
movie2gif(Mdead,{Mdead([1:2 end]).cdata},'TumourMDead.gif', ...
          'delaytime',0.1,'loopcount',0);
end
