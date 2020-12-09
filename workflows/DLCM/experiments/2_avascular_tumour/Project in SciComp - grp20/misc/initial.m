% Investigate initial conditions 

% Cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[V,R] = mesh2dual(P,E,T,'voronoi');

% Initial population
IC = 6; % Choose initial condition (1,2,3,4,5,6)
r = sqrt(P(1,:).^2+P(2,:).^2);
R1 = 0.35; % Radius of whole initial tumour
R2 = 0.1; % Radius of inner initial setup (doubly occupied, dead etc.)

if IC == 1 % circular blob of cells
    ii = find(r < 0.05); % radius of the initial blob
    U = fsparse(ii(:),1,1,[Nvoxels^2 1]);
    
elseif IC == 2 % doubly occupied voxels in center
    r1 = find(r < R1); % radius whole initial tumour
    r2 = find(r < R2); % radius of doubly occupied cells
    r1 = setdiff(r1,r2);
    U = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); 2*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]); 
    
elseif IC == 3 % dead cells in center
    r1 = find(r < R1); % radius of whole initial tumour
    r2 = find(r < R2); % radius of doubly occupied cells
    r1 = setdiff(r1,r2);
    U = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); -1*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]); 
            
elseif IC == 4 % random cell 
    r1 = find(r < R1);
    r2 = find(r < R2);
    r1 = setdiff(r1,r2);
    
    U = fsparse([r1(:); r2(:)],1, ...
                [randi([1,2],size(r1,2),1); randi([-1,0],size(r2,2),1)], ...
                [Nvoxels^2 1]); 
    
elseif IC == 5 % cell with dead center and doubly occupied voxels on boundary 
    r1 = find(r < R1);
    r2 = find(r < R2);
    r3 = find(r < R1 + 0.008);
    r1 = setdiff(r1,r2);
    r3 = setdiff(r3, union(r1,r2));
    
    U = fsparse([r1(:); r2(:); r3(:)],1, ...
                [ones(size(r1,2),1); -1*ones(size(r2,2),1); 2*ones(size(r3,2),1)], ...
                [Nvoxels^2 1]); 
            
else % "random" from previous run
    U = load('U_random.mat','U').U;
   
end

% Patch image of current cells
figure(10), clf,
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
drawnow;
