%Chrishani Jayaweera 2020-11-04
%Euler forward 



%Proliferation
U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;

u_birth__= sum(U_new<0);
if  u_birth__>0
    u_birth__
end


%Death
U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
U_new(ind_cutoff) = 0;

 u_death__= sum(U_new<0);
 if  u_death__>0
	 u_death__	
 end

% Degradation
U_deadnew(ddof) = U_deadnew(ddof) - degrade_conc*dt;
U_deadnew(U_deadnew < cutoff_deg) = 0; % remove cells below cutoff_deg

u_deg__= sum(U_new<0);
if  u_deg__>0
    u_deg__
end

% sdof_m
U_new(Adof) = U_new(Adof) + rates_sdof.*dt;

u_sdof__= sum(U_new<0);
if u_sdof__>0
    u_sdof__
end


% bdof_m
U_new(Adof) = U_new(Adof) + rates_bdof*dt;

u_bdof__= sum(U_new<0);
if u_bdof__>0
    u_bdof__
end



