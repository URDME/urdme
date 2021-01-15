
doffigure = 0;  %visualizes the dead cells in every voxel 
normalfigure = 0;   %visualizes the concentration in every voxel 
deadfigure = 0; %visualizes the sdof's and bdof's as well as the concentration in other voxels
oxyfigure = 1;
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
%%
if normalfigure==1
% create a GIF animation
Mnormal = struct('cdata',{},'colormap',{});
figure(11), clf,

Umat=full(cell2mat(Usave));
colorbar
caxis([0 max(max(Umat))])
colorlabel('Concentration of cells, U')
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = Umat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');     

    ii = find(Usave{i} == 0 & Udsave{i} > 0);
    p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    legend(p_dead,'dead')
    
    %title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    title(sprintf('Time = %d',tspan(i)));
    drawnow;
    Mnormal(i) = getframe(gcf);
end
% saves the GIF
movie2gif(Mnormal,{Mnormal([1:2 end]).cdata},'TumourMnormal.gif', ...
          'delaytime',0.1,'loopcount',0);
end

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

% %% Plot the number of cells through time
% figure(5), clf
% U_tot = sum(full(cell2mat(Usave)));
% Ud_tot = sum(full(cell2mat(Udsave)));
% plot(tspan,U_tot,'r')
% hold on
% plot(tspan,Ud_tot,'k')
% xlabel('time')
% ylabel('sum of cells in all voxels')
% legend('alive','dead')
% 
% %% Plot the number of sources through time
% figure(6), clf
% sources = sum(cell2mat(Usave) > 1);
% nonsources = sum(cell2mat(Usave) <= 1);
% plot(tspan,sources,'r')
% hold on
% plot(tspan,nonsources,'b')
% xlabel('time')
% ylabel('number of voxels')
% legend('sources','non-sources')
% 
% return;
% %%
% % investigate the time evolution of the different cell numbers
% figure(4), clf
% spsum  = @(U)(full(sum(abs(U))));
% deadsum = @(U)(full(sum(U == -1)));
% normsum = @(U)(full(sum(U == 1)));
% prolsum = @(U)(full(sum(U == 2)));
% z = cellfun(deadsum,Usave);
% w = cellfun(prolsum,Usave);
% x = cellfun(normsum,Usave);
% y = cellfun(spsum,Usave);
% p1 = plot(tspan,y);
% hold on
% p2 = plot(tspan,z,'k');
% p3 = plot(tspan,w);
% p4 = plot(tspan,x);
% p3.Color = graphics_color('vermillion');
% p4.Color = graphics_color('bluish green');
% ylim([0 max(y)]);
% xlabel('time')
% ylabel('N cells')
% legend('total', 'dead','double','single');
% 
% %%%%%%%%%%%%%%%Extra plots
% 
% %% Plot the maxium radius through time
% figure(6), clf
% plot(tspan,max_radius);
% xlabel('time')
% ylabel('max radius')
% grid on;
% 
% %% Plot the rates through time
% figure(7), clf
% rate_names = fieldnames(Ne);
% inspect_rates_norm = inspect_rates./sum(inspect_rates,1);
% bar(inspect_rates_norm','stacked','LineStyle','none') %'DisplayName',rate_names{kk});
% grid on;
% title('Relative and normalized rates')
% xlabel('time')
% ylabel('rates')
% % ticks = 
% set(gca, 'XTick', linspace(1,length(tspan),7))
% set(gca, 'XTickLabel', round(linspace(1,tspan(end),7)))
% ylim([0 1.5]);
% legend(rate_names);
% 
% 
% %% Plot Pressure
% figure(8), clf,
% Pr_ = full(U); Pr_(adof) = Pr(adof_);
% [x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,Nvoxels));
% Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
% surf(x_Pr_,y_Pr_,Pr_reshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',Pr_reshape,...
%     'EdgeColor','none');
% map_start = graphics_color('bluish green');
% map_stop = graphics_color('vermillion');
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% colormap(mymap)
% % colorbar;
% freezeColors;
% hold on;
% Pr_(adof) = 0;
% Pr_(idof) = Pr(idof_);
% Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
% surf(x_Pr_,y_Pr_,Pr_reshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',Pr_reshape,...
%     'EdgeColor','none');
% hold off;
% title('Pressure in adof(green/orange) and idof(blue)')
% map_start = [0,0,0];
% map_stop = [0,0,1];
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% caxis([-0.5 0]);
% colormap(mymap)
% %% Save the important data in a struct
% if doSave
%     saveData = struct('U', {U}, 'Usave', {Usave}, 'tspan', {tspan}, ...
%         'R', {R}, 'V', {V}, 'BC1', {BC1}, 'BC2', {BC2}, ...
%         'max_radius', {max_radius}, 'Ne', {Ne}, ...
%         'inspect_rates', {inspect_rates}, 'alpha', {alpha}, 'Pr', {Pr}, ...
%         'Adof', {Adof}, 'Nvoxels',{Nvoxels});
%     filename_saveData = "saveData/saveData_T" + Tend + ...
%         "_" + strjoin(string(fix(clock)),'-') + ".mat";
%     save(filename_saveData, 'saveData');
% end
% 
% return;
% 
% %% %% Plot U
% figure(9), clf,
% U_plt = full(U); U_plt(adof) = U(adof);
% [x_U_plt,y_U_plt] = meshgrid(linspace(-1,1,Nvoxels));
% U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
% surf(x_U_plt,y_U_plt,U_pltreshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',U_pltreshape,...
%     'EdgeColor','none');
% map_start = graphics_color('bluish green');
% map_stop = graphics_color('vermillion');
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% colormap(mymap)
% % colorbar;
% freezeColors;
% hold on;
% %U_plt(adof) = 0;
% U_plt(idof) = U(idof_);
% U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
% surf(x_U_plt,y_U_plt,U_pltreshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',U_pltreshape,...
%     'EdgeColor','none');
% hold off;
% title('Pressure in adof(green/orange) and idof(blue)')
% map_start = [0,0,0];
% map_stop = [0,0,1];
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% caxis([-0.5 0]);
% colormap(mymap)
% 
% %% %% Plot Oxygen
% figure(10), clf,
% Oxy_ = full(U); Oxy_(adof) = Oxy(adof);
% [x_Oxy_,y_Oxy_] = meshgrid(linspace(-1,1,Nvoxels));
% Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
% surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',Oxy_reshape,...
%     'EdgeColor','none');
% map_start = graphics_color('bluish green');
% map_stop = graphics_color('vermillion');
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% colormap(mymap)
% % colorbar;
% freezeColors;
% hold on;
% Oxy_(adof) = 0;
% Oxy_(idof) = Oxy(idof_);
% Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
% surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
%     'FaceAlpha','flat',...
%     'AlphaDataMapping','scaled',...
%     'AlphaData',Oxy_reshape,...
%     'EdgeColor','none');
% hold off;
% title('Pressure in adof(green/orange) and idof(blue)')
% map_start = [0,0,0];
% map_stop = [0,0,1];
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% caxis([-0.5 0]);
% colormap(mymap)
