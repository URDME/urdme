
% Report back and save time series 
if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    
    % save relevant values 
    Usave(i+1:iend) = {U};
    Udsave(i+1:iend) = {U_dead};
    
    Oxysave(i+1:iend) = {Oxy};

    bdofsave(i+1:iend) = {bdof_m};
    sdofsave(i+1:iend) = {sdof_m};
    sdofbsave{i+1:iend} = {sdof_b};

    i = iend;

    % monitor the maximum outlier cell:
    max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));
end
