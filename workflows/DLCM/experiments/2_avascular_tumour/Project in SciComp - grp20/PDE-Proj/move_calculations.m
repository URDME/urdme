
%Movement calculations
D = 1;

%sdof_m
rates_sdof = zeros(length(Adof),1);


% Pr_diff_ble = 0; 
% for ind=1:length(sdof_m_)
%     ix = sdof_m(ind);
%     ix_ = sdof_m_(ind); 
% 
%     jx_ = find(N(ix,Adof)); 
%     %Pr_diff = max(Pr(ix_)-Pr(jx_),0);
%     Pr_diff = max(Pr(ix_)-Pr(jx_),0)*(U(ix)-1);    %proportionellt mot over-occupancy
%     Pr_diff_ble= Pr_diff_ble + sum(Pr_diff);
%     
%     rates_sdof(jx_) = rates_sdof(jx_) + D*Pr_diff;
%     rates_sdof(ix_) = rates_sdof(ix_) - sum(D*Pr_diff);
% 
% end

rates_sdof__ = zeros(length(Adof),1);


[ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
%keep = find(U(Adof(jj_)) < 1);   % ...to move to
%ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
% remove any possibly remaining negative rates
Pr_diff__ = max(Pr(sdof_m_(ii))-Pr(jj_),0).*(U(sdof_m(ii))-1);    %proportionellt mot over-occupancy
grad_sdof = fsparse(ii,1,Pr_diff__*D, numel(sdof_m)); 
%moves = full(gradquotient*grad);
rates_sdof__(sdof_m_) = -gradquotient*grad_sdof;

grad_N = fsparse(jj_,1, Pr_diff__*D, numel(Adof)); 
rates_sdof__ = rates_sdof__ + gradquotient*grad_N;

rates_sdof = rates_sdof__;

% %Euler for sdof_m
% U_new(Adof) = U_new(Adof) + rates_sdof.*dt; 


% bdof_m
rates_bdof = zeros(length(Adof),1);

%check if boundary is updated. If no bdofs exist, skip bdof calculation
if (sum(bdof_m) + sum(sdof_b)) == 0
    return
elseif sum(bdof_m) == 0
    updLU = true;
    return          %!Check this expression if integrated into larger code
else
    %bdof
    %rates_bdof = zeros(length(Adof),1);

    for ind=1:length(bdof_m_)
        % movement of a boundary (singly occupied) voxel
        ix = bdof_m(ind);
        ix_ = bdof_m_(ind);

        jx_ = find(N(ix,Adof));
        % (will only move into an empty voxel:)
        jx_ = jx_(U(Adof(jx_)) == 0 & U_dead(Adof(jx_)) == 0);
        Pr_diff = max(Pr(ix_)-Pr(jx_),0);
        rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff;
        rates_bdof(ix_) = -sum(D*Pr_diff); 

    end

%     %Euler for bdof_m
%     U_new(Adof) = U_new(Adof) + rates_bdof*dt; 

    updLU = true;
end
