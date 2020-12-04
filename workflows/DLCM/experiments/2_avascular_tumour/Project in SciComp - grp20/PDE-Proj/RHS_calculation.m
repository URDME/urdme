%Pressure and oxygen calculation

if updLU
        % pressure Laplacian
        La.X = L(Adof,Adof);
        Lai = fsparse(idof_,idof_,1,size(La.X)); %remove emtpy voxels touching occupied ones
        La.X = La.X-Lai*La.X+Lai;
        [La.L,La.U,La.p,La.q,La.R] = lu(La.X,'vector');
        
        updLU = false; % assume we can reuse
end

% RHS source term proportional to the over-occupancy and BCs
Pr = full(fsparse(sdof_,1,(U(sdof)-1)./dM(sdof), ...     % Take U_dead into consideration?
    [size(La.X,1) 1]));     % RHS first...
Pr(La.q) = La.U\(La.L\(La.R(:,La.p)\Pr)); % ..then the solution

% RHS source term proportional to the over-occupancy and BCs
Oxy = full(fsparse([extdof; adof],1, ...
    [ones(size(extdof)); ...
    -cons*full(U(adof)./dM(adof))], ... 
    [size(OLa.X,1) 1]));
Oxy(OLa.q) = OLa.U\(OLa.L\(OLa.R(:,OLa.p)\Oxy));
