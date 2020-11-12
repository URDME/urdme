% Function that allocates the boundary matrix Mgamma.
% Output: Mgamma the allocated boundary matrix
% Input: p the point matrix, t the connectivity matrix 
function Mgamma = assemble_Mgamma(P,T,idof1,idof1_,LaX)
    
    % Rows and columns of T where boundary points exists (triangles on
    % the boundary idof1)
    [row,col] = find(ismember(T(1:3,:),idof1));

    % indices to unique values in T columns
    [~, ind] = unique(col, 'rows');
    % get duplicate indices i.e. the column induces of those that 
    % include two points from idof1
    duplicate_ind = setdiff(1:size(col, 1), ind);
    
    % get columns that 
    cip = zeros(2,length(duplicate_ind));
    for i = 1:length(duplicate_ind)
        cip(1,i) = T(row(duplicate_ind(i)-1),col(duplicate_ind(i)-1));
        cip(2,i) = T(row(duplicate_ind(i)),col(duplicate_ind(i)));
    end
    cip = sort(cip,1);
    [~,inds] = unique(cip(1,:));
    cip = cip(:,inds);
    
    % allocate Mgamma matrix
    Mgamma = sparse(size(LaX,1),size(LaX,2));
    
    for K = 1:size(cip,2)
          loc = cip(:,K);
          norm_val = norm(P(:,loc(1)) - P(:,loc(2)));
          MK = [2,1;1,2]/6*norm_val;
          [ind,~] = find(idof1 == loc');
          Mgamma(ind,ind) = Mgamma(ind,ind)+ MK; % add element masses to M
    end
end