      
    %proliferation-----------------------------
    ind_prol = find((Oxy > cutoff_prol));
    prol_conc = r_prol*U(ind_prol);
    
    %U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;
    
    %death--------------------------------------
    ind_die = find(Oxy < cutoff_die); %index for dying cells
    dead_conc = r_die*U(ind_die);
   
    %U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
    %U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
    
    %ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
    %U_new(ind_cutoff) = 0;
        
    %degradation--------------------------------
    degrade_conc = U_deadnew(ddof)*r_degrade; 
    %U_deadnew(ddof) = U_deadnew(ddof) - degrade_conc*dt;
    % remove degraded cells
    %U_deadnew(U_deadnew < cutoff_deg) = 0; % remove cells below cutoff_deg