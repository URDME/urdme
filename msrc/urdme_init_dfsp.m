% urdme_init_dfsp(). Creates the lookuptable (sparse matrix) for use 
%            with the dfsp solver. 
%

function fem  = urdme_init_dfsp(fem,varargin)
    dt_set=0;
    M_set=0;
    use_cache=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(fem,'comsol')
        fem  = comsol2urdme(fem.comsol);
    end
    %if ~isfield(fem,'urdme')
    %    %fem=fem2rdme(fem);
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M  = 5;
    dt = 1e-1;
    Error_Tolerance = 1e-03; %default value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    verbose=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % first check for verboseness
    if(size(varargin,1)>0 && iscell(varargin{1}))
        i=1;
        while(i<size(varargin{1},2))
            if(strcmpi(varargin{1}{i},'verbose'))
                verbose =varargin{1}{i+1};
                %fprintf('urdme_init_dfsp(): verbose=%i\n',verbose);
                break;
            end
            i=i+2;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(size(varargin,1)>0 && iscell(varargin{1}))
        i=1;
        while(i<size(varargin{1},2))
            %fprintf('%i < %i\n',i,size(varargin{1},2));
            if(strcmp(varargin{1}{i},'Error_Tolerance'))
                Error_Tolerance=varargin{1}{i+1};
                if(verbose),fprintf('\tError_Tolerance=%g\n',Error_Tolerance);end
                i=i+2;
            elseif(strcmp(varargin{1}{i},'DFSP_cache'))
                if(verbose),fprintf('checking for a cached DFSP lookup table\n');end
                dfsp_cache_filename=varargin{1}{i+1};
                if(verbose),fprintf('\t%s\n',dfsp_cache_filename);end
                use_cache=1;
                [rv,fem] = urdme_init_DFSP_read_cache(fem,dfsp_cache_filename,Error_Tolerance);
                if(rv==1),return,end
                i=i+2;
            elseif(strcmp(varargin{1}{i},'tau'))
                dt_set=1;
                dt = varargin{1}{i+1};
                if(verbose),fprintf('\ttau_d = %g\n',dt);end
                i=i+2;
            elseif(strcmp(varargin{1}{i},'max_jump'))
                M_set=1;
                M = varargin{1}{i+1};
                if(verbose),fprintf('\tMAX = %g\n',M);end
                i=i+2;
            else
                %fprintf('DFSP option not recognized:\n\t%s\n',char(varargin{1}{i}));
                i=i+2;
            end
        end
    end
    %return;
    
    %dt_set=1; M_set =1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [nspec,Ncells] = size(fem.u0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(dt_set==0 || M_set==0)
        sumDiag = -sum(diag(fem.D));
        maxDiag = max(abs(diag(fem.D)));
        if(dt_set==0 && M_set==0)
            M = urdme_init_DFSP_findM(fem,Ncells,Error_Tolerance,sumDiag,maxDiag,0);
        elseif(M_set==0)
            M = urdme_init_DFSP_findM(fem,Ncells,Error_Tolerance,sumDiag,maxDiag,dt);
        end
        if(dt_set==0)
            dt = urdme_init_DFSP_dt(fem,M,Error_Tolerance,sumDiag,maxDiag,Ncells);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if(verbose),fprintf('DFSP: using parameters tau_D=%g MAX=%g\n',dt,M);end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [DT,err] = urdme_init_DFSP_uniformization(fem,dt,M,Error_Tolerance);
%    [DT,err] = urdme_init_DFSP_ode(fem,dt,M);
   
    if(verbose),fprintf('DFSP: finished state-space exploration, error (max/avg/std) = %e/%e/%e\n',max(err),mean(err),std(err));end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fem.err = err;
    fem.DT  = DT;
    fem.Error_Tolerance = Error_Tolerance;
    % Solver options. 
    sopts(1)=dt;
    sopts(2)=M;
    fem.sopts = sopts;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(use_cache)
        if(verbose),fprintf('writing DFSP lookup table to cache');end
        urdme_init_DFSP_write_cache(fem,dfsp_cache_filename);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DT,err] =  urdme_init_DFSP_sse(fem,dt,M)
    D = fem.D;
    [Ndofs,foo] = size(fem.D);

    % 10 is an arbitrary number. Fix this! 
    % should be a function of M and the dimentionality of the system
    DT = spalloc(Ndofs,Ndofs,10*nnz(D));

    err = zeros(1,Ndofs);

    if(verbose),fprintf('Starting State-Space exploration\n');end
    sse_timer = tic;
    last_percent_reported=0;
    for i=1:Ndofs
        if(floor(i/Ndofs*100) > last_percent_reported + 4 )
            %fprintf('%g%% complete\n',floor(i/Ndofs*100));
            last_percent_reported = floor(i/Ndofs*100);
            elapsed_time = toc(sse_timer);
            time_remain = floor((1-(i/Ndofs))/((i/Ndofs)/elapsed_time));
            if(verbose>1),fprintf('%g%% complete\t\t\telapsed: %g\teta: %is\n',last_percent_reported,elapsed_time,time_remain);end
        end
       
       tier1 = i; 
       total = tier1;
       
       % Expand by reachability M steps. 
       for j=1:M
           
           tir = [];
           for k=1:numel(tier1)
              [ir,jc,s]=find(D(:,tier1(k)));       
              tir = [tir,ir'];
           end
           
          tier1 = setdiff(tir,total);
          total = union(total,tir);

       end
       tier1 = total;
       
       % tier1 now contains all states that can be touched in M steps. States
       % are sorted in ascending order in dof number. 
       
       % Assemble small FSP matrix.
       N    = numel(tier1)+1;
       DFSP = zeros(N,N); 
       DFSP(1:N-1,1:N-1) = D(tier1,tier1); 
       
       % get diagonal
       d    = diag(DFSP);
       DFSP = DFSP-diag(d);
       d    = sum(DFSP);
       
       % Add absorbing state. States in the M:th tier can go to the 
       % absorbing state. 
       
       % Difference between diagonal element in DFSP and D
       diagfull = diag(D(tier1,tier1));
       abso     = diagfull+d(1:N-1)';
       DFSP(N,1:N-1)=-abso;
       d = sum(DFSP);
       
       % Assemble diagonal into small matrix.
       DFSP = DFSP-diag(d);
       DFSP(N,N)=0;
       
       % Compute probabilities
       
       % initital condition. 
       ind = tier1 == i;
       p0  = zeros(N,1);
       p0(ind) =1;
       
       
       pdv = expm(dt*DFSP)*p0;
      % tol = 0.01*Error_Tolerance;
      % [t,p] = ode45(@rhs,[0 dt],p0,{'Abstol', tol},DFSP);
      % pdv = p(end,:);
       
       % Assemble probabilities into the table, spread the error
       % evenly to the reachable states.
       
       DT(tier1,i) = pdv(1:N-1)+pdv(N)/(N-1); 
       err(i)=pdv(N);
           
    end


end

function [DT,err] =  urdme_init_DFSP_uniformization(fem,tau,M,Error_Tolerance)

    D = fem.D;
    [Ndofs,foo] = size(fem.D);

    % 10 is an arbitrary number. Fix this! 
    % should be a function of M and the dimentionality of the system
    DT = spalloc(Ndofs,Ndofs,100*nnz(D));

    err = zeros(1,Ndofs);
     
    % Drop probabilities with lower probability than droptol relative
    % to the most likely state. 
    droptol = 0.1*Error_Tolerance;
    
    if(verbose),fprintf('Starting State-Space exploration\n');end
    sse_timer = tic;
    last_percent_reported=0;
   
    p0 = speye(Ndofs,Ndofs);
    % Uniformization
    lambda_max = max(abs(diag(D)));
    I = speye(Ndofs,Ndofs);
    A = I+D/lambda_max;
    nmax = 100;
    % tolp should be <= The DFSP tolerance \epsilon.  
    tolp=droptol;
    temp = p0;
    pdv = poisspdf(0,lambda_max*tau)*p0;
    
    totp=1;
    for i=1:nmax
        
       pp  = poisspdf(i,lambda_max*tau);
       temp = A*temp; 
       pdv = pdv+pp*temp;
       totp=totp-pp;
    
       if(totp<tolp)
           break;
       end
    end
    
    % Truncate 
    for i=1:Ndofs
        
        if(floor(i/Ndofs*100) > last_percent_reported + 4 )
            %fprintf('%g%% complete\n',floor(i/Ndofs*100));
            last_percent_reported = floor(i/Ndofs*100);
            elapsed_time = toc(sse_timer);
            time_remain = floor((1-(i/Ndofs))/((i/Ndofs)/elapsed_time));
            if(verbose>1),fprintf('%g%% complete\t\t\telapsed: %g\teta: %is\n',last_percent_reported,elapsed_time,time_remain);end
        end
       
       maxp = max(pdv(:,i)); 
       [tier1,jj,ss] = find(pdv(:,i)/maxp>droptol);  
       ss = pdv(tier1,i);
       % Assemble probabilities into the DFSP table, spread the error
       % evenly to the reachable states.
       ee = 1-norm(ss,1);
       err(i)=ee;
       DT(tier1,i) = ss/norm(ss,1); 
       
           
    end


end

function [DT,err] =  urdme_init_DFSP_ode(fem,tau)

    D = fem.D;
    [Ndofs,foo] = size(fem.D);

    % 10 is an arbitrary number. Fix this! 
    % should be a function of M and the dimentionality of the system
    DT = spalloc(Ndofs,Ndofs,100*nnz(D));

    err = zeros(1,Ndofs);
     
    % Drop probabilities with lower probability than droptol relative
    % to the most likely state. 
    droptol = 1e-4;
    
    if(verbose),fprintf('Starting State-Space exploration\n');end
    sse_timer = tic;
    last_percent_reported=0;
   
    p0 = speye(Ndofs,Ndofs);
    I = speye(Ndofs,Ndofs);
    dt = tau;
    A1 = I+0.5*dt*D';
    A2 = I-0.5*dt*D';
        
    % Truncate 
    for i=1:Ndofs
        
        if(floor(i/Ndofs*100) > last_percent_reported + 4 )
            %fprintf('%g%% complete\n',floor(i/Ndofs*100));
            last_percent_reported = floor(i/Ndofs*100);
            elapsed_time = toc(sse_timer);
            time_remain = floor((1-(i/Ndofs))/((i/Ndofs)/elapsed_time));
            if(verbose>1),fprintf('%g%% complete\t\t\telapsed: %g\teta: %is\n',last_percent_reported,elapsed_time,time_remain);end
        end
        
       % pdv = A1\(A2*p0(:,i)); 
       %[t,pdv] = ode45(@ode_rhs,[0 tau],p0(:,i),[],D);

        
       
       maxp = max(pdv); 
       [tier1,jj,ss] = find(pdv/maxp>droptol);  
       ss = pdv(tier1,i);
       % Assemble probabilities into the DFSP table, spread the error
       % evenly to the reachable states.
       ee = 1-norm(ss,1);
       err(i)=ee;
       DT(tier1,i) = ss/norm(ss,1); 
       
       i    
    end


end

function y = ode_rhs(t,x,D)

    y = D*x;    

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = urdme_init_DFSP_findM(fem,Ncells,Error_Tolerance,sumDiag,maxDiag,dt_in)
    % DFSP is faster than NSM if the expected number of diffusional 
    %  events (a_diff*tau_D) is greater than the number of voxels (K).
    %            a_diff * tau_D > K
    %            sumDiag * tau_D > K
    %            tau_D > K / sumDiag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt =  max([ Ncells/sumDiag dt_in]);
    if(verbose),fprintf('urdme_init_DFSP_findM: dt = %g\n',dt);end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maximum_M = 25;
    minum_M = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find optimum MAX parameter
    M=-1;
    for i= minum_M:maximum_M
        tmp_err = 1-sum(poisspdf(0:i,maxDiag*dt));
        if(verbose),fprintf('\tM=%g E[err]=%g  (maxDiag=%g)\n',i,tmp_err,full(maxDiag));end
        if(tmp_err < Error_Tolerance)
            M = i;
            %fprintf('MAX=%g\n',M);
            %error('stop');
            return;
        end
    end
    if(M<0), error('Could not find MAX parameter'),end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find tau_d
function dt = urdme_init_DFSP_dt(fem,M,Error_Tolerance,sumDiag,maxDiag,Ncells)
    dt_low=0;
    dt_low_set=0;
    %
    dt_high = Inf;
    dt_high_set=0;
    %
    dt = Ncells/sumDiag;  % initial guess
    dt_set = 0;
    %
    dt_max  = min(diff(fem.tspan));
    %
    while(dt_set ~=1)
        if ~isfinite(dt)
            error('non-finite dt');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_err = 1-sum(poisspdf(0:M,maxDiag*dt));
        if(verbose),fprintf('DFSP tau_D=%g\t\terror=%e\t\ttol=%e\n',dt,tmp_err,Error_Tolerance);end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if( tmp_err < Error_Tolerance )
            %fprintf('\ttoo small\n');
            dt_low_set=1;
            dt_low=dt;
            if(dt_high_set)
                dt_new = (dt_high-dt_low)/2 + dt_low;
                if(dt_new < 1.1*dt)
                    dt=dt_new;
                    break;
                else
                    dt=dt_new;
                end
            else
                dt = dt*10;
            end
            if dt > dt_max
                dt=dt_max;
                break;
            end
        else
            %fprintf('\ttoo big\n');
            dt_high_set=1;
            dt_high=dt;
            if(dt_low_set)
                dt = (dt_high-dt_low)/2;
            else
                dt = dt/10;
            end
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rv,fem] = urdme_init_DFSP_read_cache(fem,dfsp_cache_filename,Error_Tolerance)
    rv=0;
    fid=fopen(dfsp_cache_filename); %does the file exist?
    if(fid==-1)
        return;
    else
        fclose(fid);
        verified=0;
        %%%%%
        tmp_fem = load(dfsp_cache_filename);
        if(isfield(tmp_fem,'urdme'))
            tmp_u = tmp_fem.urdme;
            if(...
            isfield(tmp_u,'Error_Tolerance') &&...
            isfield(tmp_u,'DT') &&...
            isfield(tmp_u,'sopts') &&...
            isfield(tmp_u,'err') ...
            )
                if((tmp_fem.urdme.Error_Tolerance == Error_Tolerance) &...
                   (tmp_fem.mesh.p == fem.mesh.p) )
                    verified=1;
                end
            end
        end
        %%%%%
        if verified==1
            fem.urdme.DT = tmp_fem.urdme.DT;
            fem.urdme.err = tmp_fem.urdme.err;
            fem.urdme.sopts = tmp_fem.urdme.sopts;
            fem.urdme.Error_Tolerance = tmp_fem.urdme.Error_Tolerance;
            rv=1;
            return
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function urdme_init_DFSP_write_cache(fem,dfsp_cache_filename)
    if(verbose),fprintf('saving DFSP matrix as %s\n',dfsp_cache_filename);end
    save(dfsp_cache_filename,'-STRUCT','fem');
end



end %function urdme_init_DFSP

