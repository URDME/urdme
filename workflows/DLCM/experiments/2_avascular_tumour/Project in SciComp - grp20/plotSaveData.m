%% Plot all saved data and compare
% Johannes Dufva 2020-11-06

% Get all saved .mat files from saveData-folder
DirList = dir(fullfile('saveData/', '*.mat'));
Data = cell(1, length(DirList));
figrows = ceil(sqrt(length(DirList)));

%% Plot population appearance

% figure(1);

for k = 1:length(DirList)
    load(['saveData/' DirList(k).name]);
    loadSaveData;
    t_show = find(tspan >= 150, 1); % length(Usave);
    if(tspan(end) >= 150)
%         subplot(figrows, ceil(length(DirList)/figrows),k);
        figure(k)
        patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
                    'EdgeColor','none');
        hold on,
        axis([-1 1 -1 1]); axis square, axis off
        ii = find(Usave{t_show} == 1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('bluish green'));
        ii = find(Usave{t_show} == 2);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',graphics_color('vermillion'));
        ii = find(Usave{t_show} == -1);
        patch('Faces',R(ii,:),'Vertices',V, ...
            'FaceColor',[0 0 0]);
        title(sprintf('Time = %d, Ncells = %d,\n BC1 = %d, BC2 = %d', ...
            tspan(t_show),full(sum(abs(Usave{t_show}))), BC1, BC2));
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