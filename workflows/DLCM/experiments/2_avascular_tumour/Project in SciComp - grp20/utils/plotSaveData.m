%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
folder = 'saveData/2020-11-30, basePr=1_T500_IC5/';
DirList = dir(fullfile(folder, '*.mat'));
Data = cell(1, length(DirList));
figrows = ceil(sqrt(length(DirList)));

%% Plot population appearance
%close all;
for k = 1:length(DirList)
    load([folder, DirList(k).name]);
    t_show = length(tspan);
    if tspan(end) == 500 || tspan(end) == 0
        figure('Name',"CellPlot_" + DirList(k).name(10:end-4));
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
                    'EdgeColor','none');
        hold on,
%         axis([-1 1 -1 1]); axis square%, axis off
        axis([-1 1 -1 1]*0.55); axis square, %axis off
        ii = find(Usave{t_show} == 1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('bluish green'));
        ii = find(Usave{t_show} == 2);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('vermillion'));
        ii = find(Usave{t_show} == -1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
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
    if tspan(end) == 600 || tspan(end) == 0
        figure;
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
        ylim([0 max(y)]);
        title(sprintf('alpha = %d',alpha));
        xlabel('time')
        ylabel('N cells')
        legend('total', 'dead','double','single');
    end
end

%% Plot rates
for k = 1:length(DirList)
    load([folder DirList(k).name]);
    if tspan(end) == 10 || tspan(end) == 0
        figure('Name',"RatePlot_" + DirList(k).name(10:end-4));
        plotRates;
    end
end

%% Plot pressure
% close all;
for k = 1:length(DirList)
%     subplot(figrows, ceil(length(DirList)/figrows),k);  
    load([folder DirList(k).name]);
    if tspan(end) == 500 || tspan(end) == 0
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