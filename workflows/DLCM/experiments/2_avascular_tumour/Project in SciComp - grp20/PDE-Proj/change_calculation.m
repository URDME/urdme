% Change calculation of proliferation, death and degradation

%proliferation-----------------------------
ind_prol = find((Oxy > cutoff_prol));
prol_conc = r_prol*U(ind_prol);

%death--------------------------------------
ind_die = find(Oxy < cutoff_die); %index for dying cells
dead_conc = r_die*U(ind_die);

%degradation--------------------------------
degrade_conc = U_deadnew(ddof)*r_degrade; 
