% Patch image of current bdof_m
figure(11), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
patch('Faces',R(bdof_m,:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'));
drawnow;