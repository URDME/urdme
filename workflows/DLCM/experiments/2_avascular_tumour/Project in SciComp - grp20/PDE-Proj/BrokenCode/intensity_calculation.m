
% Intensities of possible events
moveb = rates_bdof(bdof_m_);
moves = rates_sdof(sdof_m_);

%New calculation
%proliferation
birth = (ind_prol>0)*r_prol; 

%death
death = (ind_die>0)*r_die; 

%degradation
degrade = (r_degrade*(U_deadnew(ddof)>0));

%Total intensity 
%intens = [birth; death; degrade; moves; moveb];






