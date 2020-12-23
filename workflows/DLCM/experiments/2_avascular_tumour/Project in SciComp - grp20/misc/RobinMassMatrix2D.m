function R = RobinMassMatrix2D(p,e)
    np = size(p,2); % number of nodes
    ne = size(e,2); % number of boundary edges
    R = sparse(np,np); % allocate boundary matrix
    
    for E = 1:ne
        loc2glb = e(1:2,E); % boundary nodes
        x = p(1,loc2glb); % node x-coordinates
        y = p(2,loc2glb); % node ylen
        len = sqrt((x(1)-x(2))^2+(y(1)-y(2))^2); % edge length
        RE = [2 1; 1 2]*len/6; % edge boundary matrix
        R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + RE;
    end
end