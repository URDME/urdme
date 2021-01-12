% Patch image of current cells
figure(666);
hold on,
axis([-1 1 -1 1]); axis square, axis off
ii = find(U <= 1 & U > 0);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'));
ii = find(U > 1);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',graphics_color('vermillion'));
ii = find(U < 0);
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 0]);
ii = idof1;
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 0 1]);
ii = idof2;
patch('Faces',R(ii,:),'Vertices',V, ...
    'FaceColor',[0 1 1]);
title(sprintf('Time = %d',tt));
drawnow;