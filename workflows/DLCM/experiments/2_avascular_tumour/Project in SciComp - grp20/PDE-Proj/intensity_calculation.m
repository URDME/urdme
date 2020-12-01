    % intensities of possible events
    %   % (1) moving boundary DOFs
      [ii,jj_] = find(N(bdof_m,Adof)); % neighbours...
      keep = find(U(Adof(jj_)) == 0);  % ...to move to
      ikeep= ii(keep);
      ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1); %???varför??
      % remove any possibly remaining negative rates
      grad = fsparse(ii,1,max(Pr(bdof_m_(ii))-Pr(jj_),0).* ...
                     1,...%Drate_(2*VU(Adof(jj_))+1), ... % (U(Adof(jj_)) = 0)
                     numel(bdof_m));
      moveb = full(gradquotient*grad);
    %
    % (2) also certain sources may move by the same physics
    [ii,jj_] = find(N(sdof_m,Adof)); % neighbours...
    keep = find(U(Adof(jj_)) < 1);   % ...to move to
    ii = reshape(ii(keep),[],1); jj_ = reshape(jj_(keep),[],1);
    % remove any possibly remaining negative rates
    grad = fsparse(ii,1,max(Pr(sdof_m_(ii))-Pr(jj_),0).* ...
        1, ... %Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1)
        numel(sdof_m)); % (abs as U could be -1)
    moves = full(gradquotient*grad);
    
    %     moves=moveb;
    % (3) proliferation/death/degradation rates
%     r_prol*(U(Adof) > 0 & U(Adof) < 2).*(Oxy(Adof) > cutoff_prol);
    birth = full(r_prol*(U(Adof) > 0 & U(Adof) < 2).*(Oxy(Adof) > cutoff_prol)); %U(Adof) < 2 sätter gräns för när den får proliferate, 1?
    total_birth = sum(birth);
    birth = total_birth/total_birth * birth;
    birth(isnan(birth)) = 0;
    % (as we get some 0/0 terms if total_birth == 0);
    
    death = full(r_die*(U(Adof) > 0).*(Oxy(Adof) < cutoff_die));
    degrade = full(r_degrade*(U_dead > 0));
    %
    %   %Gillespies algorithm
    %    intens = [moveb; moves; birth; abs(death); abs(degrade)];
    intens = [birth; death; degrade; moves; moveb];
    
    