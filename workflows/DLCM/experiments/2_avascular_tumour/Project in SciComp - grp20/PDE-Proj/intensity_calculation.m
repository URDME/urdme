
% bdof_m
moveb = rates_bdof(bdof_m_);

% sdof_m
moves = rates_sdof(sdof_m_);

% birth
birth = (ind_prol>0)*r_prol; 

% death
death = (ind_die>0)*r_die; 

% degradation
degrade = (r_degrade*(U_deadnew(ddof)>0));

% total
intens = [birth; death; degrade; moves; moveb];

%%%#########
intens_print__ = [sum(birth); sum(death); sum(degrade); sum(moves); sum(moveb)]; 





