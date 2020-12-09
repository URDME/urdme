function [P,E,T,gradquotient] = weird_mesh(Nvoxels)
% WEIRD_MESH

    gd = [3 4 -1 1 1 -1 -1 -1 1 1]';
    sf = 'SQ1';
    ns = char(sf)';
    G = decsg(gd,sf,ns);
    [P,E,T] = poimeshSUPER(G,Nvoxels-1);

    % figure;
    % sznt = size(T,2);
    % for nt = 1:sznt
    %     x = P(1,T(1:3,nt));
    %     y = P(2,T(1:3,nt));
    %     patch(x,y,nt/sznt);
    %     text(mean(x),mean(y),num2str(nt));
    % end
    % colorbar

    if mod(Nvoxels,2) == 0
        ind_vec = [2*ones(1, floor((Nvoxels-1)/2)^2 + ceil((Nvoxels-1)/2)^2), 1, ...
            2*ones(1, -1 + floor((Nvoxels-1)/2)^2 + ceil((Nvoxels-1)/2)^2)];
        start = -1;
        k = 0;
    else
        ind_vec = [2*ones(1,(Nvoxels-1)/2), 1, 2*ones(1,(Nvoxels-1)/2)];
        ind_vec(1) = [];
        ind_vec(end) = 3;
        k = 1;
        start = 0;
    end

    temp = [];
    T_2 = T;
    while start < size(T,2)
            start = start + ind_vec(1+mod(k,length(ind_vec)));
            temp = [temp,start];
            k = k + 1;
    end

    split_ind = size(T,2)/2;
    temp_1 = temp(temp <= split_ind);
    temp_2 = temp(temp > split_ind);
    T_2(3, temp_1) = T_2(3,temp_1) - 1;
    T_2(1, temp_2) = T_2(1,temp_2) + 1;
    
    T = T_2;
    % distance between two voxel midpoints
    h = 2/(Nvoxels-1);
    D1 = h;

    % edge length
    L1 = h;

    % quotient
    gradquotient = L1/D1; % (= 1 for Cartesian mesh)
end
