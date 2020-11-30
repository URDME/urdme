% Patch image of current cells
% figure;
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
ii = find(U == 1);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'));
ii = find(U == 2);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('vermillion'));
ii = find(U == -1);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0]);
title(sprintf('Time = %d, Ncells = %d',tt,full(sum(abs(U)))));
drawnow;
% P(1,idof2)
% P(2,idof2)