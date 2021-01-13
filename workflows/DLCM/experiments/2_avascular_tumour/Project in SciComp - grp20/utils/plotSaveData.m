%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
folder = 'saveData/2021-01-08, FEMfinalSetup_T1000_IC1/';
% folder = 'saveData/';
DirList = dir(fullfile(folder, '*.mat'));
Data = cell(1, length(DirList));
figrows = ceil(sqrt(length(DirList)));

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
    t_show = length(tspan);
    if tspan(end) == 2000 || tspan(end) == 0
        figure('Name',"CellPlot_" + DirList(k).name(10:end-4));
% %         patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
% %                     'EdgeColor','none');
        hold on,
        U_show = Usave{t_show};
%         axis([-1 1 -1 1]); axis square%, axis off
        axis([-1 1 -1 1]*0.75); axis square, axis off
        ii = find(U_show == 1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('bluish green'));
        ii = find(U_show == 2);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('vermillion'));
        ii = find(U_show == -1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
%         ii = idof(ismember(idof,find(~VU)));
%         patch('Faces',R(ii,:),'Vertices',V, ...
%             'FaceColor',[0 0 1]);
%         ii = idof(ismember(idof,find(VU)));
%         patch('Faces',R(ii,:),'Vertices',V, ...
%             'FaceColor',[1 1 0]);
%         title(sprintf('Time = %d, Ncells = %d \n' , ...
%             tspan(t_show),full(sum(abs(Usave{t_show})))));
        title(sprintf('Time = %d, Ncells = %d \n \\alpha = %e' , ...
            tspan(t_show),full(sum(abs(Usave{t_show}))), alpha));
%         title(sprintf('Time = %d, Ncells = %d \n Ne.moveb2 = %d' , ...
%             tspan(t_show),full(sum(abs(Usave{t_show}))), Ne.moveb2));
        hold off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow;
    end
end

%% Plot population appearance (over time)
alphaToShow = 1e-2;
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
    t_show_vec = 1:20:length(tspan);
    if alpha == alphaToShow
        for t_show = t_show_vec
            fig = figure('Name',tspan(t_show) + "_" + DirList(k).name(10:end-22));
            U_show = Usave{t_show};%"CellPlotTIME_t + " 
            hold on
%             axis([-1 1 -1 1]); axis square, axis off
            axis([-1 1 -1 1]*0.75); axis square, axis off
            ii = find(U_show == 0);
            patch('Faces',R(ii,:),'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
            'EdgeColor','none');
            ii = find(U_show == 1);
            patch('Faces',R(ii,:),'Vertices',V, ...
                'FaceColor',graphics_color('bluish green'));
            ii = find(U_show == 2);
            patch('Faces',R(ii,:),'Vertices',V, ...
                'FaceColor',graphics_color('vermillion'));
            ii = find(U_show == -1);
            patch('Faces',R(ii,:),'Vertices',V, ...
                'FaceColor',[0 0 0]);
%             title(sprintf('Time = %d\n', ...
%                 tspan(t_show)));
            drawnow;
%             text(-0.60,0.60,sprintf('[t = %d]', tspan(t_show)),...
%                 'Fontsize',14,'FontName','Gill Alt One MT Italic');
            hold off;
            set(gca,'LooseInset',get(gca,'TightInset'));
        end
    end
end

%% Plot number of cell plots

spsum  = @(U)(full(sum(abs(U))));
deadsum = @(U)(full(sum(U == -1)));
normsum = @(U)(full(sum(U == 1)));
prolsum = @(U)(full(sum(U == 2)));

for k = 1:length(DirList)
    load([folder DirList(k).name]);
    if tspan(end) == 1000 || tspan(end) == 230
        figure('Name',"EvolutionPlot_" + DirList(k).name(10:end-4));
        z = cellfun(deadsum,Usave);
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

figure('Name',"TOTALPlot_" + folder(10:end-1));
hold on;
tToShow = Tend;
alphaToShow = [1e-4, 1e-1];
plotStyle = {'-o','-s','-^','-d'};
ymax = 1;
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    if tspan(end) >= tToShow && ismember(alpha,alphaToShow)
        y = cellfun(spsum,Usave);%mod(k-1,length(plotStyle))
        p1 = plot(tspan,y,plotStyle{1+(alpha == 1e-4)},...
            'MarkerIndices',1:10:length(y),...
            'Displayname',sprintf('\\alpha = %1.0e', alpha), ...
            'Linewidth', 1.3);
        if max(y) > ymax
            ymax = max(y);
        end
    end
end
plot(tspan,ones(length(y),1)*y(1),'Displayname','Initial cells', ...
            'Linewidth', 1.0,'LineStyle','--','Color','k');
ymin = min(y);
% title(sprintf('Total cells at t = %d',tToShow));
xlabel('time','Fontsize',12)
xticks([0,200,400,600,Tend])
ylabel('N cells','Fontsize',12)
% yticks([1400:300:3200])
ax = gca;
ax.FontSize = 12;
legend('Fontsize',10);
% ylim([ymin ymax]);
grid on;
% hold off;

%% Plot number of total cells together 2 MISC
spsum  = @(U)(full(sum(abs(U))));

figure('Name',"TOTALPlot_" + folder(10:end-1));
hold on;
tToShow = Tend;
alphaToShow = [1e-4, 1e-1, 1e+4];
plotStyle = {'-o','-s'};
plotColor = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980]};
ymax = 1;
for k = 1:length(DirList)
    fname = DirList(sort_k(k)).name;
    load([folder fname]);
    if tspan(end) >= tToShow && ismember(alpha,alphaToShow)
        y = cellfun(spsum,Usave);
        p1 = plot(tspan,y,plotStyle{1+(string(fname(end-6:end-4)) == "FEM")},...
            'Color', ((50+k)/50)*plotColor{1+(string(fname(end-6:end-4)) == "FEM")},...
            'MarkerIndices',1:10:length(y),...
            'Displayname',sprintf([fname(end-6:end-4) ' \\alpha = %1.0e'], alpha), ...
            'Linewidth', 1.3);
        if max(y) > ymax
            ymax = max(y);
        end
    end
end
% plot(tspan,ones(length(y),1)*y(1),'Displayname','Initial cells', ...
%             'Linewidth', 1.0,'LineStyle','--','Color','k');
ymin = min(y);
% title(sprintf('Total cells at t = %d',tToShow));
xlabel('time','Fontsize',12)
xticks([0,Tend])
ylabel('N cells','Fontsize',12)
% yticks([1400:300:3200])
ax = gca;
ax.FontSize = 12;
legend('Fontsize',10);
plots=get(gca, 'Children');
legend(plots([6,4,2,5,3,1]));
% ylim([ymin ymax]);
grid on;
% hold off;

%% Plot Ne
load([folder DirList(1).name]);
fnames = fieldnames(Ne);
alphaToShow = [1e-4, 1e-2, 1e-1, 1e+4];
Ne_vector = zeros(length(fnames),length(alphaToShow));
kk = 1;
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    if tspan(end) >= 400 && ismember(alpha,alphaToShow)
        Ne_vector(:,kk) = cell2mat(struct2cell(Ne))';
        kk = kk + 1;
    end
end
% Remove moveb2 elements (due to its relatively huge impact)
% Ne_vector(2,:) = [];
% fnames(2) = [];

% Plot stacked bar plot
figure('Name',"NePlot_" + folder(10:end-4));
b = bar(Ne_vector','LineStyle','none');

% Plot settings
grid on;
title('Ne without moveb2');
xlabel('alpha')
ylabel('rates')
xticks(1:length(alphaToShow))
xticklabels(alphaToShow);
legend(fnames);
%% Plot rates
for k = 1:length(DirList)
    load([folder DirList(k).name]);
    if tspan(end) == 2000 || tspan(end) == 0
        figure('Name',"RatePlot_" + DirList(k).name(10:end-4));
        plotRates;
    end
end

%% Plot pressure
% close all;
for k = 1:length(DirList)
%     subplot(figrows, ceil(length(DirList)/figrows),k);  
    load([folder DirList(k).name]);
    if tspan(end) == 1000 || tspan(end) == 0
        figure('Name',"PrPlot_" + DirList(k).name(10:end-4));
        plotPressure;
    end
end

%% Print all open figures
% Observe that you need to name all the figures that you want to save
openFigures = findobj('Type', 'figure');
for k = 1:length(openFigures)
    figNumb = openFigures(k).Number;
    figHandle = "-f" + figNumb;
    figName = openFigures(k).Name;
    print(figHandle, '-r300', "images/" + figName,'-painters','-depsc'); %    '-painters'
end