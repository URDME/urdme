if i==1
    Oxysave{1} = Oxy(Adof);
end

% report back håll some koll
if tspan(i+1) < tt+dt
    iend = i+find(tspan(i+1:end) < tt+dt,1,'last');
    Usave(i+1:iend) = {U};
    Udsave(i+1:iend) = {U_dead};
    Idofsave(i+1:iend) = {U(Idof)};

    Oxysave(i+1:iend) = {Oxy(Adof)};
    bdofsave(i+1:iend) = {bdof_m};
    sdofsave(i+1:iend) = {sdof_m};

    i = iend;

    % monitor the maximum outlier cell:
    max_radius = sqrt(max(P(1,adof).^2+P(2,adof).^2));

    % the number of cells
    num_cells = sum(abs(U));
end
