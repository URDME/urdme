%shitForMakingInitialPic
%%
% Some test

fig420 = figure(420);
fig420.Name = 'alt_StartStruct_idof3_shifted';
ah = gca;
hold on;

pdemesh(P,E,T);
plot(ah, P(1,:), P(2,:),'.', 'MarkerSize', 10, 'Color', 'b');

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

plt_alpha = 0.5;
% ii = setdiff(1:size(R,1),idof);
% patch(ah,'Faces',R(ii,:),'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
%     'EdgeColor','none');
% axis([-1 1 -1 1 0 1]), axis off,  axis tight
axis([-1 1 -1 1 0 1]*0.7), axis off, axis square
ii = find(Pr == 1);
patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'),'FaceAlpha',plt_alpha);%,'EdgeColor','none');
ii = find(Pr == 2);
patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('vermillion'),'FaceAlpha',plt_alpha);%,'EdgeColor','none');
ii = find(Pr == -1);
patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0],'FaceAlpha',plt_alpha);%,'EdgeColor','none');
ii = idof1;
patch(ah,'Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0,0,1],'FaceAlpha',plt_alpha-0.1);%,'EdgeColor','none');
% title('Start out structure')
drawnow;
% view(3);

for k = 1:length(Adof)
    text(ah,P(1,Adof(k)),P(2,Adof(k)),sprintf('p_{%d}',k),'Color','w','FontSize',11,'HorizontalAlignment','center');
end


% Making plot tight and keep colors right
outerpos = ah.OuterPosition;
ti = ah.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ah.Position = [left bottom ax_width ax_height];
fig420.InvertHardcopy = 'off';