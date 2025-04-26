% Function that sets initial conditions for tumour
function U_initial = set_init(IC,R1,R2,P,Nvoxels,r_p)
% IC: Initial condition (1,2,3,4,5,6)
%   IC1 --> orginally IC, single occupied cells
%   IC2 --> doubly occupied voxels in center
%   IC3 --> dead cells in center
%   IC4 --> random; dead/empty in center surrounded by single/doubly occupied
%   IC5 --> tumour with dead center and doubly occupied voxels on boundary
%   IC6 --> tumour with empty center and --||-- 
%   IC7 --> tumour with dead center, rest is singly occupied
%   IC8 --> tissue outide circular tumour, sources on ring outside tumour
% P: point matrix from mesh
% R1: Radius of whole initial tumour
% R2: Radius of inner initial setup (doubly occupied, dead etc.)
% r_p: additional geometry parameter. Its function depends on which IC:
%   IC 1: blob radius
%   IC 5,6: thickness, in number of cells, of proliferating ring


r = sqrt(P(1,:).^2+P(2,:).^2);
h = 2/(Nvoxels-1);

if IC == 1 % circular blob of cells
    ii = find(r < r_p); % radius of the initial blob %0.05
    U_initial = fsparse(ii(:),1,1,[Nvoxels^2 1]);
    
elseif IC == 2 % doubly occupied voxels in center
    r1 = find(r < R1); % radius whole initial tumour
    r2 = find(r < R2); % radius of doubly occupied cells
    r1 = setdiff(r1,r2);
    U_initial = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); 2*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]); 
    
elseif IC == 3 % dead cells in center
    r1 = find(r < R1); % radius of whole initial tumour
    r2 = find(r < R2); % radius of dead cells
    r1 = setdiff(r1,r2);
    U_initial = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); -1*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]); 
    
elseif IC == 4 % random cell 
    r1 = find(r < R1);
    r2 = find(r < R2);
    r1 = setdiff(r1,r2);
    
    U_initial = fsparse([r1(:); r2(:)],1, ...
                [randi([1,2],size(r1,2),1); randi([-1,0],size(r2,2),1)], ...
                [Nvoxels^2 1]); 
            
elseif IC == 5 % cell with dead center and doubly occupied voxels on boundary 
    r1 = find(r < R1);
    r2 = find(r < R2);
    r3 = find(r < R1 + r_p*h);
    r1 = setdiff(r1,r2);
    r3 = setdiff(r3, union(r1,r2));
    
    U_initial = fsparse([r1(:); r2(:); r3(:)],1, ...
                [ones(size(r1,2),1); -1*ones(size(r2,2),1); 2*ones(size(r3,2),1)], ...
                [Nvoxels^2 1]);
            
elseif IC == 6 % cell with empty center and doubly occupied voxels on boundary 
    r1 = find(r < R1);
    r2 = find(r < R2);
    r3 = find(r < R1 + r_p*h);
    r1 = setdiff(r1,r2);
    r3 = setdiff(r3, union(r1,r2));
    
    % Unecessary to initialize zeros in sparse matrix?
    U_initial = fsparse([r1(:); r2(:); r3(:)],1, ...
                [ones(size(r1,2),1); zeros(size(r2,2),1); 2*ones(size(r3,2),1)], ...
                [Nvoxels^2 1]);
elseif IC == 7 % singly occupied voxels only, dead center
    r1 = find(r < R1);
    r2 = find(r < R2);
    r3 = find(r < R1 + h);
    r1 = setdiff(r1,r2);
    r3 = setdiff(r3, union(r1,r2));
    
    U_initial = fsparse([r1(:); r2(:); r3(:)],1, ...
                [ones(size(r1,2),1); -1*ones(size(r2,2),1); 1*ones(size(r3,2),1)], ...
                [Nvoxels^2 1]);
elseif IC == 8 % tissue outside circular tumour, sources at R2 ring
    r1 = find(r > R1);
    r2 = find(r > R2);
    r3 = find(r > R2 + h);
    r2 = setdiff(r2, r3);
    d = 0.05;
    r4 = find( P(1,:) > -d &  P(1,:) < d & P(2,:) > R2 & P(2,:) < R2 + d);
    r2 = union(r2,r4);
    
    U_initial = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); 2*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]);
elseif IC == 9 % tissue outside circular tumour, random sources outside R2
    r1 = find(r > R1);
    r2 = find(r > R2);
    r1 = setdiff(r1,r2);
    
    U_initial = fsparse([r1(:); r2(:)],1, ...
                [ones(size(r1,2),1); 2*ones(size(r2,2),1)], ...
                [Nvoxels^2 1]);
elseif IC == 10 %  square tumour
    r1 = find(P(1,:) < R1 & P(1,:) > -R1);
    r2 = find(P(2,:) < R2 & P(2,:) > -R2);
    r1 = intersect(r1,r2);
    
    U_initial = fsparse([r1(:)],1, ...
                [ones(size(r1,2),1)], ...
                [Nvoxels^2 1]);
            
elseif IC == 11
    
else
    error('Initial Condition Type not valid');
    
end