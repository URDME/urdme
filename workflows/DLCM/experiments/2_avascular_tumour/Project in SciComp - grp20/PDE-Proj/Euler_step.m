%Chrishani Jayaweera 2020-11-04
%Euler forward 

%Proliferation
U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;

%Death
U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
U_new(ind_cutoff) = 0;

% Degradation
U_deadnew(ddof) = U_deadnew(ddof) - degrade_conc*dt;
U_deadnew(U_deadnew < cutoff_deg) = 0; % remove cells below cutoff_deg

% sdof_m
U_new(Adof) = U_new(Adof) + rates_sdof.*dt;
% bdof_m
U_new(Adof) = U_new(Adof) + rates_bdof*dt;