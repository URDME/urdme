
    % classify the DOFs
    adof = find(U|U_dead); % all filled voxels 
    
    % singularly occupied voxels on the boundary: 
    bdof_m = find(N*(U ~= 0 | U_dead ~= 0) < neigh & (U > cutoff_bdof & ...
        U <= 1));
    
    sdof = find(U > 1); % voxels with 2 cells
    sdof_b = find(N*(U~=0 | U_dead ~=0) < neigh & (U > 1));
    
    % voxels with more than concentration 1 in them which may move, 
    % with a voxel containing less number of cells next to it:

    %sdof_m = intersect(find(sum(N.*(U'+U_dead')<(U+U_dead)& logical(N),2)), find(U > 1));
    %sdof_m = find(sum(N.*(U'+U_dead')<(U+U_dead)&(N&N),2).*(U>1));
    sdof_m = find(sum(N.*(U')<(U)&(N&N),2).*(U>1));

    %FIX!!!! Look at total cell-concetration in neighbour (living and dead)
    
    Idof = (N*(U ~= 0) > 0 & U == 0); % empty voxels touching occupied ones
    idof1 = find(Idof & ~VU);         % "external" OBC1
    idof2 = find(Idof & VU);          % "internal" OBC2
    idof = find(Idof);
    ddof = find(U_dead > 0);   %degrading voxels
    
    % "All DOFs" = adof + idof, like the "hull of adof"
    Adof = [adof; idof];
    % The above will be enumerated within U, a Nvoxels^2-by-1 sparse
    % matrix. Determine also a local enumeration, eg. [1 2 3
    % ... numel(Adof)].

    Adof_ = (1:numel(Adof))';
    [bdof_m_,sdof_,sdof_m_,idof1_,idof2_,idof_,adof_, sdof_b_,ddof_] = ...
        map(Adof_,Adof,bdof_m,sdof,sdof_m,idof1,idof2,idof,adof,sdof_b,ddof);
    