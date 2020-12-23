%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
folder = 'saveData/2020-12-16, scaleL_neigh_NEWaLai_idof3_IC5/';
% folder = 'saveData/';
DirList = dir(fullfile(folder, '*.mat'));
Data = cell(1, length(DirList));
figrows = ceil(sqrt(length(DirList)));

% Get sorted alpha values
a = 0;
alphas = [];
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
        alphas = [alphas, alpha];
end
sorted_alphas = sort(alphas);
[~,sort_k] = ismember(sorted_alphas,alphas);

%% Plot population appearance
%close all;
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
    t_show = length(tspan);
    if tspan(end) == 400 || tspan(end) == 0
        figure('Name',"CellPlot_" + DirList(k).name(10:end-4));
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
                    'EdgeColor','none');
        hold on,
        U_show = Usave{t_show};
        axis([-1 1 -1 1]); axis square%, axis off
%         axis([-1 1 -1 1]*0.65); axis square, %axis off
        ii = find(U_show == 1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('bluish green'));
        ii = find(U_show == 2);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('vermillion'));
        ii = find(U_show == -1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
        ii = idof(ismember(idof,find(~VU)));
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 1]);
        ii = idof(ismember(idof,find(VU)));
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[1 1 0]);
        title(sprintf('Time = %d, Ncells = %d \n alpha = %d' , ...
            tspan(t_show),full(sum(abs(Usave{t_show}))), alpha));
%         title(sprintf('Time = %d, Ncells = %d \n Ne.moveb2 = %d' , ...
%             tspan(t_show),full(sum(abs(Usave{t_show}))), Ne.moveb2));

        drawnow;
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

figure('Name',"TOTALPlot_");
hold on;
tToShow = 400;
ymax = 1;
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    if tspan(end) == tToShow
        y = cellfun(spsum,Usave);
        p1 = plot(tspan,y,'Displayname',sprintf('\\alpha = %d', alpha), ...
            'Linewidth', 1.0);
        if max(y) > ymax
            ymax = max(y);
        end
    end
end
ymin = min(y);
title(sprintf('Total cells at t = %d',tToShow));
xlabel('time')
ylabel('N cells')
legend();
% ylim([ymin ymax]);
grid on;
hold off;

%% Plot Ne
load([folder DirList(1).name]);
fnames = fieldnames(Ne);
Ne_vector = zeros(length(fnames),length(DirList));
for k = 1:length(DirList)
    load([folder DirList(sort_k(k)).name]);
    if tspan(end) == 400 || tspan(end) == 0
        Ne_vector(:,k) = cell2mat(struct2cell(Ne))';
        
    end
end
% Remove moveb2 elements (due to its relatively huge impact)
Ne_vector(2,:) = [];
fnames(2) = [];

% Plot stacked bar plot
figure('Name',"NePlot_" + folder);
b = bar(Ne_vector','stacked','LineStyle','none');

% Plot settings
grid on;
title('Ne without moveb2');
xlabel('alpha')
ylabel('rates')
xticks(1:length(sorted_alphas))
xticklabels(sorted_alphas);
legend(fnames);
%% Plot rates
for k = 1:length(DirList)
    load([folder DirList(k).name]);
    if tspan(end) == 400 || tspan(end) == 0
        figure('Name',"RatePlot_" + DirList(k).name(10:end-4));
        plotRates;
    end
end

%% Plot pressure
% close all;
for k = 1:length(DirList)
%     subplot(figrows, ceil(length(DirList)/figrows),k);  
    load([folder DirList(k).name]);
    if tspan(end) == 400 || tspan(end) == 0
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
    print(figHandle, '-r300', "images/" + figName, '-dpng');
end