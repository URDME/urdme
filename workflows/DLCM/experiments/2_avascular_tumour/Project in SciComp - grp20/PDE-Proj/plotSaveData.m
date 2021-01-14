%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
folder = 'testShit/2020-01-12, T=150_exp=0/';
% folder = 'testShit/';
DirList = dir(fullfile(folder, '*.mat'));

% Get sorted alpha values
a = 0;
alphas = [];
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
        alphas = [alphas, alpha*k];
end
sorted_alphas = sort(alphas);
[~,sort_k] = ismember(sorted_alphas,alphas);

%% Plot population appearance
%close all;
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
    t_show = find(tspan == 150);
        figure('Name',"CellPlot_" + DirList(k).name(10:end-4));
        set(gca,'color','none');
        Umat=full(cell2mat(Usave));
%         colorbar('southoutside')
        caxis([0 max(max(Usave{t_show}))])
%         colorlabel('Concentration of cells, U')
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off
        ii = find(Usave{t_show}>0);
        c = full(Usave{t_show});
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceVertexCData',c(ii),'FaceColor','flat');     
        ii = find(Usave{t_show} == 0 & Udsave{t_show} > 0);% ,'Edgecolor', 'none'
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
        hold off;
%         title(sprintf('Time = %d, Ncells = %d \n \\alpha = %1.0e' , ...
%             tspan(t_show),full(sum(Usave{t_show})) + full(sum(Udsave{t_show})), alpha));       
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow;
end

%% Plot population appearance (over time)
alphaToShow = 1e-4;
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
    t_show_vec = 1:25:length(tspan);
    if alpha == alphaToShow
        for t_show = t_show_vec
            figure('Name',"CellPlotTIME_t" + t_show + "_" + DirList(k).name(10:end-4));
            set(gca,'color','none');
            Umat=full(cell2mat(Usave));
            colorbar('southoutside')
            caxis([0 max(max(Usave{t_show}))])
            colorlabel('Concentration of cells, U')
    %         patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    %             'EdgeColor','none');
            hold on,
            axis([-1 1 -1 1]); axis square, axis off
            ii = find(Usave{t_show}>0);
            c = full(Usave{t_show});
            patch('Faces',R(ii,:),'Vertices',V, ...
                'FaceVertexCData',c(ii),'FaceColor','flat');%,'Edgecolor', 'none'
            ii = find(Usave{t_show} == 0 & Udsave{t_show} > 0);
            patch('Faces',R(ii,:),'Vertices',V, ...
                'FaceColor',[0 0 0]);
            hold off;
            title(sprintf('Time = %d, Ncells = %d \n \\alpha = %1.0e' , ...
                tspan(t_show),full(sum(Usave{t_show})) + full(sum(Udsave{t_show})), alpha));       
            set(gca,'LooseInset',get(gca,'TightInset'));
            drawnow;
        end
    end
end

%% Plot number of cell plots

spsum  = @(U)(full(sum(abs(U))));
deadsum = @(U)(full(sum(abs(U))));
normsum = @(U)(full(sum(U <= 1 & U > 0)));
prolsum = @(U)(full(sum(U > 1)));

for k = 1:length(DirList)
    load([folder DirList(k).name]);
    if tspan(end) == 200
        figure('Name',"EvolutionPlot_" + DirList(k).name(10:end-4));
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
%         ylim([0 max(y)]);
%         ylim([0 1500]);
        xlim([0 150]);
        title(sprintf('alpha = %d',alpha));
        xlabel('time')
        ylabel('N cells')
        legend('total', 'dead','double','single');
        grid on;
        hold off;
    end
end

%% Plot number of total cells together
spsum  = @(U)(full(sum(abs(U))));

figure('Name',"TOTALPlot_" + folder(10:end-4));
hold on;
Tend = 100;
tToShow = find(tspan == Tend);
alphaToShow = [1e-4, 1e+4, Inf];
plotStyle = {'-o','-s','-^','-d'};
ymax = 1;
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    if ismember(alpha,alphaToShow)
        y = cellfun(spsum,Usave(1:tToShow));
        y = y + cellfun(spsum,Udsave(1:tToShow));
        p1 = plot(tspan(1:tToShow),y,plotStyle{1+mod(k-1,length(plotStyle))},...
            'MarkerIndices',1:10:length(y),...
            'Displayname',sprintf('\\alpha = %1.0e', alpha), ...
            'Linewidth', 1.3);
        if max(y) > ymax
            ymax = max(y);
        end
    end
end
% plot(tspan(1:tToShow),ones(length(y),1)*y(1),'Displayname','Initial cells', ...
%             'Linewidth', 1.0,'LineStyle','--','Color','k');
ymin = min(y);
% title(sprintf('Total cells at t = %d',tToShow));
xlabel('time','Fontsize',12)
xticks([0,150])
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
ylabel('N cells','Fontsize',12)
% yticks([1400:300:3200])
ax = gca;
ax.FontSize = 12;
legend('Fontsize',12);
% ylim([ymin ymax]);
grid on;
hold off;

%% Plot rates
for k = 1:length(DirList)
    load([folder DirList(k).name]);
    figure('Name',"RatePlot_" + DirList(k).name(10:end-4));
    plotRates2;
end
%% Plot pressure
% close all;
for k = 1:length(DirList)
    load([folder DirList(k).name]);
    figure('Name',"PrPlot_" + DirList(k).name(10:end-4));
    plotPressure;
end

%% Plot pressure diff in bars
ymax=0;
fig1 = figure('Name',"PressureDiffPlot_" + DirList(1).name(10:end-4));
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    subplot(2,1,k);
%     fig1.Position = [0 300 500 400];
    [iii,jjj_] = find(N(idof,adof)); % neighbours...
    grad = max(Pr(jjj_) - Pr(idof_(iii)),0);
    bar(grad);
    % Fix data tip labels to show pressure between which nodes
    labels = compose('%d -> %d', [idof(iii),jjj_]);
    dcm_obj = datacursormode(fig1);
    set(dcm_obj,'UpdateFcn',{@datacursor,labels})
    % set ylim and title
    xlabel('Connections from tumour to boundary','Fontsize',12);
    ylabel('\nabla P','Fontsize',12);
    xticks([0:100:400])
    if ymax < max(grad)
        ymax = max(grad);
    end
    if k < 2
      set(gca,'xticklabel',{[]})
      xlabel('');
    end   
    ylim([0 ymax]);
    grid on;
%     title(sprintf('sum bars = %d \n \\alpha = %1.0e', sum(grad),alpha));
    title(sprintf('\\alpha = %1.0e',alpha),'Fontsize',12);
    ax = gca;
    ax.FontSize = 11;
end
%% Print all open figures
% Observe that you need to name all the figures that you want to save
openFigures = findobj('Type', 'figure');
for k = 1:length(openFigures)
    figNumb = openFigures(k).Number;
    figHandle = "-f" + figNumb;
    figName = openFigures(k).Name;
    print(figHandle, '-r300', "testshit/imgs/" + figName,'-painters','-depsc'); %    '-painters'
end