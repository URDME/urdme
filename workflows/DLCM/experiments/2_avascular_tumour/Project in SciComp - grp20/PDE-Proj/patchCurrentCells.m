% Patch image of current cells
figure(67);
set(gca,'color','none');
Umat=full(cell2mat(Usave));
colorbar('southoutside')
caxis([0 max(max(U))])
colorlabel('Concentration of cells, U')
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off

ii = find(U>0);
c = full(U);
patch('Faces',R(ii,:),'Vertices',V,'FaceVertexCData',c(ii),'FaceColor','flat');     

ii = find(U == 0 & U_dead > 0);
p_dead = patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0]);
title(sprintf('Time = %d',tt));
drawnow;