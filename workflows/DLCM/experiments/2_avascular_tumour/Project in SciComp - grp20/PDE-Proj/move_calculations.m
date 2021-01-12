
%sdof_m
rates_sdof = zeros(length(Adof),1);

for ind=1:length(sdof_m_)
    ix = sdof_m(ind);
    ix_ = sdof_m_(ind); 

    jx_ = find(N(ix,Adof)); 
    Pr_diff = max(Pr(ix_)-Pr(jx_),0);%*(U(ix)-1);    %proportionellt mot over-occupancy
    
    rates_sdof(jx_) = rates_sdof(jx_) + D*Pr_diff;
    rates_sdof(ix_) = rates_sdof(ix_) - sum(D*Pr_diff);
end

% bdof_m
rates_bdof = zeros(length(Adof),1);

%check if boundary is updated. If no bdofs exist, skip bdof calculation
if (sum(bdof_m) + sum(sdof_b)) ~= 0
    for ind=1:length(bdof_m_)
        % movement of a boundary (singly occupied) voxel
        ix = bdof_m(ind);
        ix_ = bdof_m_(ind);

        jx_ = find(N(ix,Adof));
        % (will only move into an empty voxel:)
        jx_ = jx_(U(Adof(jx_)) == 0 & U_dead(Adof(jx_)) == 0);
        Pr_diff = max(Pr(ix_)-Pr(jx_),0);%*U(ix);    %proportionellt mot over-occupancy

        rates_bdof(jx_) = rates_bdof(jx_) + D*Pr_diff;
        rates_bdof(ix_) = -sum(D*Pr_diff); 
    end
        updLU = true;
end
