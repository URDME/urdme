function Mgamma = assemble_Mgamma2(p,t)

    np = size(p,2); % Number of nodes 
    Mgamma = sparse(np,np); % Allocate boundary mass matrix 

    % Go through every node to find the edge elements    
    for N = 1:np

        % Find column number (triangle) in which node N is part of 
        [~, col] = find(t(1:3,:)==N);

        % Find other nodes in these triangles
        tri = t(1:3, col);

        % Find neighbours to N
        N_neigh = unique(tri(tri~=N));

        % Calculate edge lengths from N to its neighbours
        len = sqrt((p(1,N) - p(1,N_neigh)).^2 + ((p(2,N) - p(2,N_neigh)).^2));

        % Add contribution to boundary mass matrix
        Mgamma(N,N_neigh) = (1/6).*len;

    end
end