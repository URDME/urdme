% Function that allocates the boundary matrix Mgamma.
% Output: Mgamma the allocated boundary matrix
% Input: p the point matrix, t the connectivity matrix 
function Mgamma = assemble_Mgamma(P,T,idof1,idof1_,)
    np = size(p,2); % Number of nodes 
    nt = size(t,2); % Number of boundary edges
    Mgamma = sparse(np,np); % Allocate mass matrix 
    
    [row,col] = find(ismember(T(1:3,:),idof1));

    % indices to unique values in col
    [~, ind] = unique(col, 'rows');
    % duplicate indices
    duplicate_ind = setdiff(1:size(col, 1), ind);
    % % duplicate values
    % duplicate_value = col(duplicate_ind);

    cont_int_points = zeros(2,length(duplicate_ind));

    for i = 1:length(duplicate_ind)
        cont_int_points(1,i) = T(row(duplicate_ind(i)-1),col(duplicate_ind(i)-1));
        cont_int_points(2,i) = T(row(duplicate_ind(i)),col(duplicate_ind(i)));
    end
    cont_int_points = sort(cont_int_points,1);
    [~,inds] = unique(cont_int_points(1,:));
    cont_int_points = cont_int_points(:,inds);
    
    for K = 1:nt
%         loc2glb = t(1:3,K); 
%         x = p(1,loc2glb); % x-coordinates of triangle nodes
%         y = p(2,loc2glb); % y-coordinates of triangle nodes
%         area = polyarea(x,y); % Area of triangle 
%         MK = [2 1 1; 1 2 1; 1 1 2]/12*area; % element mass matrix 
%         Mgamma(loc2glb,loc2glb) = Mgamma(loc2glb,loc2glb)+ MK; % add element masses to M
    end
end