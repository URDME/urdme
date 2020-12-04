
% intensities of possible events
%   % (1) moving boundary DOFs
[ii,jj_] = find(N(bdof_m,Adof)); % neighbours... ??bdofs tat souldnt be bdofs, or pressure not calced correctly
keep = find(U(Adof(jj_)) == 0);  % ...to move to
ikeep= ii(keep);
ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1); 
% remove any possibly remaining negative rates
% bdof_m_(ii);
% Pr(bdof_m_(ii));
% max(Pr(bdof_m_(ii))-Pr(jj_),0).*D;
% numel(bdof_m)
grad = fsparse(ii,1,max(Pr(bdof_m_(ii))-Pr(jj_),0).*D,numel(bdof_m));
moveb = full(gradquotient*grad);
rates_bdof;

%moveb = (rates_bdof>0)*rates_bdof;

%
% (2) also certain sources may move by the same physics
[ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
keep = find(U(Adof(jj_)) < 1);   % ...to move to
ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
% remove any possibly remaining negative rates
grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0).* ...
    D, numel(sdof_m)); 
moves = full(gradquotient*grad);
rates_sdof;

%! sum(rates_sdof(sdof_m_)) = sum(moves)
%! rates_sdof(sdof_m_)) = moves
% figure out if we can rewrite the for-loop with the help of these lines
% fundera ur detta rknas ut

%moves = (rates_sdof>0)*rates_sdof;

% % (3) proliferation/death/degradation rates
% birth = full(r_prol*(U(Adof) > 0 & U(Adof) < 2).*(Oxy(Adof)...
%     > cutoff_prol)); %U(Adof) < 2 sätter gräns för när den får proliferate, 1?
% total_birth = sum(birth);
% birth = total_birth/total_birth * birth;
% birth(isnan(birth)) = 0;
% % (as we get some 0/0 terms if total_birth == 0);
% 
% death = full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die));
% degrade = full(r_degrade*(U_dead > 0));
% % 
% % 

% New calculation
%proliferation
%ind_prol = find(Oxy > cutoff_prol);
birth = (ind_prol>0)*r_prol; 

%death
%ind_die = find(Oxy < cutoff_die); %index for dying cells
death = (ind_die>0)*r_die; 

%degradation
degrade = (r_degrade*(U_deadnew(ddof)>0));

% total
intens = [birth; death; degrade; moves; moveb];

