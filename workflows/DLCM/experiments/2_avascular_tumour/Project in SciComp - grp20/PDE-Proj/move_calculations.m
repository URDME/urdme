    D = 1; 
    %sdof_m
    rates_sdof = zeros(length(Adof),1);
    rates_test = zeros(length(Adof),1);

    for ind=1:length(sdof_m_)
%         Ne.moves = Ne.moves+1;
        ix = sdof_m(ind);
        ix_ = sdof_m_(ind); %index i Adof

        jx_ = find(N(ix,Adof)); %index i Adof
        Pr_diff = max(Pr(ix_)-Pr(jx_),0);
        rates_sdof(jx_) = rates_sdof(jx_) + D*Pr_diff;
        rates_sdof(ix_) = rates_sdof(ix_) - sum(D*Pr_diff);

    end

    %Euler for sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt; %action_vec*

    %bdof
    rates_bdof = zeros(length(Adof),1);
%     action_vec =zeros(length(Adof),1);
%     Ne.moveb = Ne.moveb+1;
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

    %Euler for bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt; %action_vec*

