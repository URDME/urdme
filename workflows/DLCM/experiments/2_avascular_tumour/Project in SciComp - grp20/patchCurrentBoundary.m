% Patch image of current boundary
figure(12), clf,
patch('Faces',R,'Vertices',V,'FaceColor',[0.9 0.9 0.9], ...
    'EdgeColor','none');
hold on,
axis([-1 1 -1 1]); axis square, axis off
patch('Faces',R([idof1;idof3],:),'Vertices',V, ...
    'FaceColor',graphics_color('bluish green'));
drawnow;