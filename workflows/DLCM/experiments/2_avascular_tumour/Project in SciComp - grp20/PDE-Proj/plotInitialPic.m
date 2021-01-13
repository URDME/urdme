%shitForMakingInitialPic
%%
% Some test

fig420 = figure(420); clf
fig420.Name = 'what';
fig420.Position = [1000 300 600 500];
ah = gca;
hold on;

% pdemesh(P,E,T);
% ii = 1:size(P,2);
% plot(ah, P(1,ii), P(2,ii),'.', 'MarkerSize', 10, 'Color', 'b');

% hat_val = 0.4;
% Pr_ = full(zeros(size(Pr))); Pr_(idof1(1)) = hat_val; Pr_(idof1(7)) = hat_val;
% Pr_reshape = reshape(Pr_, Nvoxels, Nvoxels);
% [x_Pr_,y_Pr_] = meshgrid(linspace(-1,1,size(Pr_reshape,1)));
% pdemesh(P,E,T);
% surf(ah, x_Pr_,y_Pr_,Pr_reshape,'FaceAlpha',0,'Marker','.','MarkerSize',10,'LineWidth',1.0,'EdgeAlpha',0.5);
% alim([-1 1])
% map_start = graphics_color('bluish green');
% map_stop = graphics_color('vermillion');
% xx = linspace(0,1,10);
% map_matrix = map_start' + xx.*(map_stop' - map_start');
% mymap = map_matrix';
% caxis([mean(Pr_)-2*std(Pr_) mean(Pr_)+2*std(Pr_)]);
% colormap(mymap)

plt_alpha = 0.9;
% ii = setdiff(1:size(R,1),idof);
% patch(ah,'Faces',R(ii,:),'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%     'EdgeColor','none');
% axis([-1 1 -1 1 0 1]), axis square, axis off
axis([-1 1 -1 1 0 1]*0.7), axis square, axis off
ii = find(U < 1 & U > 0);
p = patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'),'FaceAlpha',plt_alpha);%,'EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ii = find(U > 1);
p = patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('vermillion'),'FaceAlpha',plt_alpha);%,'EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ii = find(U < 0);
p = patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0],'FaceAlpha',plt_alpha);%,'EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ii = idof;
p = patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0,0,1],'FaceAlpha',plt_alpha-0.1);%,'EdgeColor','none');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% title('Start out structure')
drawnow;
% view(3);

% ii = find(Pr);
% plot(ah, P(1,ii), P(2,ii),'square', 'MarkerSize', 10, 'Color', 'k','Displayname',sprintf('\\Omega'));
% ii = idof1;
% plot(ah, P(1,ii), P(2,ii),'.', 'MarkerSize', 10, 'Color', 'k','Displayname',sprintf('\\Gamma'));
% legend('FontSize', 14);


for k = 1:length(idof)
    text(ah,P(1,idof(k)),P(2,idof(k)),sprintf('p_{%d}',idof(k)),'Color','w','FontSize',11,'HorizontalAlignment','center');
end


% Making plot tight and keep colors right
% outerpos = ah.OuterPosition;
% ti = ah.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ah.Position = [left bottom ax_width ax_height];
% set(gca,'LooseInset',get(gca,'TightInset'));

% set(gcf, 'PaperSize', [6.25 6.25]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 6.25 6.25]);
% 
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [6.25 6.25]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 6.25 6.25]);
% fig420.InvertHardcopy = 'off';