    
    %proliferation------------------------------
%     Ne.birth = Ne.birth+1;
%     birth_count = birth_count+1;
    

    ind_prol = find((Oxy > cutoff_prol));
    prol_conc = r_prol*U(ind_prol);
    U_new(ind_prol)=U_new(ind_prol)+prol_conc*dt;
    
    %death--------------------------------------
    % full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die))
    ind_die = find(Oxy < cutoff_die); %index for dying cells
    
%     Oxy(ind_die);
    dead_conc = r_die*U(ind_die);
    U_new(ind_die) = U_new(ind_die) - dead_conc*dt;
    U_deadnew(ind_die) = U_deadnew(ind_die) + dead_conc*dt;
    
    ind_cutoff =  find(U_new < cutoff_remain & (Oxy < cutoff_die));%*check dead/alive instead?
    U_new(ind_cutoff) = 0;
    
    %Ne.death = Ne.death+1;
    
    %degradation--------------------------------
    %Ne.degrade = Ne.degrade+1;   
    U_deadnew(ddof) = U_deadnew(ddof)*(1 - r_degrade*dt); %tidssteg litet->kan inte bli negativt
    % remove degraded cells
    U_deadnew(U_deadnew < cutoff_deg) = 0; %ta bort för små