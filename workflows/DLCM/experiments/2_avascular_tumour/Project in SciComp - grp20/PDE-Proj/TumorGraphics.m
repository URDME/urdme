
load('AvascularTumor_workspace.mat')


%% 

% create a GIF animation
% figure(6)
% population appearance sdof, bdof visualization
Mdof = struct('cdata',{},'colormap',{});
figure(3), clf,

Umat=cell2mat(Usave);
cmat = full(Umat/max(max(Umat)));
colorbar
caxis([0 1])
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
        
%     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 0.5 & Usave{i} <= 1);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 1 & Usave{i} <= 1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 0]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 1.5 & Usave{i} <= 2);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 1]); %[1 0 Usave{i}/max(Usave{i})]
%     
%     ii = find(Usave{i} > 2 );
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 1 0]); %[1 0 Usave{i}/max(Usave{i})]
    
    patch('Faces',R(bdofsave{i},:),'Vertices',V, ...
        'FaceColor','cyan'); %[1 0 Usave{i}/max(Usave{i})] 
    
    patch('Faces',R(sdofsave{i},:),'Vertices',V, ...
        'FaceColor','magenta'); %[1 0 Usave{i}/max(Usave{i})]      
        
%     ii = find(Usave{i} > 3 );
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 0 0]); %[1 0 Usave{i}/max(Usave{i})]
    
    ii = find(Usave{i} == 0 & Udsave{i} >0);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    drawnow;
    Mdof(i) = getframe(gcf);
%     pause(3)
end
% 
% % saves the GIF
% movie2gif(Mdof,{Mdof([1:2 end]).cdata},'animations/TumourMdof.gif', ...
%           'delaytime',0.1,'loopcount',0);

%%
if normalfigure==1
% create a GIF animation

% population appearance normal
Mnormal = struct('cdata',{},'colormap',{});
figure(11), clf,

Umat=cell2mat(Usave);
cmat = full(Umat/max(max(Umat)));
colorbar
caxis([0 1])
for i = 1:numel(Usave)
    
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Usave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
    
    %     ii = find(Usave{i}>0 & Usave{i} <= 1);
    %     cvec = full(Usave{i}(ii)/max(Usave{i}(ii)));
    %     patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData', cvec);% ...
    %         %'FaceColor',[0 0 Usave{i}/max(Usave{i})]);    

    %     ii = find(Usave{i} > 0 & Usave{i} <=0.5);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 0 1]); %; [0 1 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 0.5 & Usave{i} <= 1);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 0.5 1]); %; [0 1 Usave{i}/max(Usave{i})]

%     ii = find(Usave{i} > 1);
%     cvec = full(Usave{i}(ii)/max(Usave{i}(ii)));
%     patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData', cvec);% ...
%         %'FaceColor',[0 0 Usave{i}/max(Usave{i})]);   
    
    %     ii = find(Usave{i} > 1 & Usave{i} <= 1.5);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0 1 0.5]); %; [0 1 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 1.5 & Usave{i} <= 2);
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[0.5 1 0.2]); %[1 0 Usave{i}/max(Usave{i})]
    %     
    %     ii = find(Usave{i} > 2 );
    %     patch('Faces',R(ii,:),'Vertices',V, ...
    %         'FaceColor',[1 1 0.5]); %[1 0 Usave{i}/max(Usave{i})]

    ii = find(Usave{i} == 0 & Udsave{i} >0);
    patch('Faces',R(ii,:),'Vertices',V, ...
        'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d, Nbdof = %d',tspan(i),full(sum(abs(Usave{i}))),length(bdofsave{i})));
    drawnow;
    Mnormal(i) = getframe(gcf);
%     pause(2)
end


%%
if deadfigure==1
% dead appearance
Mdead = struct('cdata',{},'colormap',{});
figure(10), clf,

Udmat=cell2mat(Udsave);
cmat = full(Udmat/max(max(Udmat)));
caxis([1 2])
colorbar;
%     set( h, 'YDir', 'reverse' );
colormap 'gray'

for i = 1:numel(Udsave)
    patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
        'EdgeColor','none');
    hold on,
    axis([-1 1 -1 1]); axis square, axis off
    
    ii = find(Udsave{i}>0);
    c = cmat(ii,i);
    patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c,'FaceColor','flat');    
    
    
%     ii = find(Udsave{i} > 0 & Udsave{i} <=0.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[1 1 1]); %; [0 1 Usave{i}/max(Usave{i})]
%     
%     ii = find(Udsave{i} > 0.5 & Udsave{i} <= 1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0.5 0.5 0.5]); %; [0 1 Usave{i}/max(Usave{i})] 
%   
%     ii = find(Udsave{i} >1.5);
%     patch('Faces',R(ii,:),'Vertices',V, ...
%         'FaceColor',[0 0 0]);%'FaceVertexCData',color,'FaceColor','flat');
    
    title(sprintf('Time = %d, Ncells = %d',tspan(i),full(sum(abs(Usave{i})))));
    drawnow;
    Mdead(i) = getframe(gcf);
end

%investigate the time evolution of the different cell numbers
figure(4), clf
spsum  = @(U)(full(sum(abs(U))));
deadsum = @(U)(full(sum(U_dead == -1)));
normsum = @(U)(full(sum(U == 1)));
prolsum = @(U)(full(sum(U == 2)));

z = cellfun(deadsum,Udsave);
w = cellfun(prolsum,Usave);
x = cellfun(normsum,Usave);
y = cellfun(spsum,Usave);
p1 = plot(tspan,y);
hold on
p2 = plot(tspan,z,'k');
p3 = plot(tspan,w);
p4 = plot(tspan,x);
p3.Color = graphics_color('vermillion');
p4.Color = graphics_color('bluish green');
ylim([0 max(y)]);
xlabel('time')
ylabel('N cells')
legend('total', 'dead','double','single');

return;


%%%%%%%%%%%%%%%Extra plots

%% Plot the maxium radius through time
figure(5), clf
plot(tspan,max_radius);
xlabel('time')
ylabel('max radius')
grid on;

%% Plot the rates through time
figure(6), clf
rate_names = fieldnames(Ne);
inspect_rates_norm = inspect_rates./sum(inspect_rates,1);
bar(inspect_rates_norm','stacked','LineStyle','none') %'DisplayName',rate_names{kk});
grid on;
title('Relative and normalized rates')
xlabel('time')
ylabel('rates')
% ticks = 
set(gca, 'XTick', linspace(1,length(tspan),7))
set(gca, 'XTickLabel', round(linspace(1,tspan(end),7)))
ylim([0 1.5]);
legend(rate_names);


%% Plot Pressure
figure(7), clf,
Pr_ = full(U); Pr_(adof) = Pr(adof_);
[x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,Nvoxels));
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
surf(x_Pr_,y_Pr_,Pr_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Pr_reshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
Pr_(adof) = 0;
Pr_(idof) = Pr(idof_);
Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
surf(x_Pr_,y_Pr_,Pr_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Pr_reshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)
%% Save the important data in a struct
if doSave
    saveData = struct('U', {U}, 'Usave', {Usave}, 'tspan', {tspan}, ...
        'R', {R}, 'V', {V}, 'BC1', {BC1}, 'BC2', {BC2}, ...
        'max_radius', {max_radius}, 'Ne', {Ne}, ...
        'inspect_rates', {inspect_rates}, 'alpha', {alpha}, 'Pr', {Pr}, ...
        'Adof', {Adof}, 'Nvoxels',{Nvoxels});
    filename_saveData = "saveData/saveData_T" + Tend + ...
        "_" + strjoin(string(fix(clock)),'-') + ".mat";
    save(filename_saveData, 'saveData');
end

return;

%% %% Plot U
figure(9), clf,
U_plt = full(U); U_plt(adof) = U(adof);
[x_U_plt,y_U_plt] = meshgrid(linspace(-1,1,Nvoxels));
U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
surf(x_U_plt,y_U_plt,U_pltreshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',U_pltreshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
%U_plt(adof) = 0;
U_plt(idof) = U(idof_);
U_pltreshape = reshape(U_plt, Nvoxels, Nvoxels);
surf(x_U_plt,y_U_plt,U_pltreshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',U_pltreshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)

%% %% Plot Oxygen
figure(10), clf,
Oxy_ = full(U); Oxy_(adof) = Oxy(adof);
[x_Oxy_,y_Oxy_] = meshgrid(linspace(-1,1,Nvoxels));
Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Oxy_reshape,...
    'EdgeColor','none');
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
colormap(mymap)
% colorbar;
freezeColors;
hold on;
Oxy_(adof) = 0;
Oxy_(idof) = Oxy(idof_);
Oxy_reshape = reshape(Oxy_, Nvoxels, Nvoxels);
surf(x_Oxy_,y_Oxy_,Oxy_reshape,...
    'FaceAlpha','flat',...
    'AlphaDataMapping','scaled',...
    'AlphaData',Oxy_reshape,...
    'EdgeColor','none');
hold off;
title('Pressure in adof(green/orange) and idof(blue)')
map_start = [0,0,0];
map_stop = [0,0,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([-0.5 0]);
colormap(mymap)


