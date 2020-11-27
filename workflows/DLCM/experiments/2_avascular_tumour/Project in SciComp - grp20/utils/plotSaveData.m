%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
DirList = dir(fullfile('saveData/2020-11-23/', '*.mat'));
Data = cell(1, length(DirList));
figrows = ceil(sqrt(length(DirList)));

%% Plot population appearance
%close all;
for k = 1:length(DirList)
    load(['saveData/2020-11-23/' DirList(k).name]);
    t_show = length(tspan);
    if tspan(end) > 10 || tspan(end) == 0
        figure('Name',"CellPlot_" + DirList(k).name(10:end-4));
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
                    'EdgeColor','none');
        hold on,
%         axis([-1 1 -1 1]); axis square, axis off
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
        title(sprintf('Time = %d, Ncells = %d', ...
            tspan(t_show),full(sum(abs(Usave{t_show})))));
        drawnow;
    end
end

%% Plot number of cell plots

spsum  = @(U)(full(sum(abs(U))));
deadsum = @(U)(full(sum(U == -1)));
normsum = @(U)(full(sum(U == 1)));
prolsum = @(U)(full(sum(U == 2)));

figure(2);
for k = 1:length(DirList)
%     subplot(figrows, ceil(length(DirList)/figrows),k);
    figure(k)
    load(['saveData/' DirList(k).name]);
    loadSaveData;
    
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
    title(sprintf('BC1 = %d, BC2 = %d', ...
        BC1, BC2));
    xlabel('time')
    ylabel('N cells')
    legend('total', 'dead','double','single');
end

%% Plot pressure
%close all;
for k = 1:length(DirList)
%     subplot(figrows, ceil(length(DirList)/figrows),k);  
    load(['saveData/' DirList(k).name]);
    if tspan(end) > 10 || tspan(end) == 0
        figure('Name',"PrPlot_" + DirList(k).name(10:end-4));
        plotPressure;
    end
end

%% Print all open figures

openFigures = findobj('Type', 'figure');
for k = 1:length(openFigures)
    figNumb = openFigures(k).Number;
    figHandle = "-f" + figNumb;
    figName = openFigures(k).Name;
    print(figHandle, '-r300', "images/" + figName, '-dpng');
end