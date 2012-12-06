% urdme_init_dfsp(). Creates the lookuptable (sparse matrix) for use 
%            with the dfsp solver. 
%

function fem  = urdme_init_dfsp(fem,varargin)
    dt_set=0;
    M_set=0;
    use_cache=0;
    if isfield(fem,'comsol')
        umod = comsol2urdme(fem.comsol);
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %if ~isfield(fem,'urdme')
    %    fem=comsol2urdme(fem);
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Error_Tolerance = 1e-02; %default value
    if(size(varargin,1)>0 && iscell(varargin{1}))
        i=1;
        while(i<size(varargin{1},2))
            %fprintf('%i < %i\n',i,size(varargin{1},2));
            if(strcmp(varargin{1}{i},'Error_Tolerance'))
                Error_Tolerance=varargin{1}{i+1};
                fprintf('\tError_Tolerance=%g\n',Error_Tolerance);
                i=i+2;
            elseif(strcmp(varargin{1}{i},'DFSP_cache'))
                fprintf('checking for a cached DFSP lookup table\n');
                dfsp_cache_filename=varargin{1}{i+1};
                fprintf('\t%s\n',dfsp_cache_filename)
                use_cache=1;
                [rv,fem] = urdme_init_DFSP_read_cache(fem,dfsp_cache_filename,Error_Tolerance);
                if(rv==1),return,end
                i=i+2;
            elseif(strcmp(varargin{1}{i},'tau'))
                dt_set=1;
                dt = varargin{1}{i+1};
                fprintf('\ttau_d = %g\n',dt);
                i=i+2;
            elseif(strcmp(varargin{1}{i},'max_jump'))
                M_set=1;
                M = varargin{1}{i+1};
                fprintf('\tMAX = %g\n',M);
                i=i+2;
            else
                %fprintf('DFSP option not recognized:\n\t%s\n',char(varargin{1}{i}));
                i=i+2;
            end
        end
    end
    %return;
    
    M  = 5;
    dt = 1;
    dt_set=1; M_set =1;
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
    fprintf('DFSP: using parameters tau_D=%g MAX=%g\n',dt,M);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [DT,err] = urdme_init_DFSP_uniformization(fem,dt,M);
    fprintf('DFSP: finished state-space exploration, error (max/avg/std) = %e/%e/%e\n',max(err),mean(err),std(err));
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
        fprintf('writing DFSP lookup table to cache');
        urdme_init_DFSP_write_cache(fem,dfsp_cache_filename);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end %function urdme_init_DFSP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [DT,err] =  urdme_init_DFSP_uniformization(fem,tau,M)

    D = fem.D;
    [Ndofs,foo] = size(fem.D);

    % 10 is an arbitrary number. Fix this! 
    % should be a function of M and the dimentionality of the system
    DT = spalloc(Ndofs,Ndofs,100*nnz(D));

    err = zeros(1,Ndofs);
     
    % Drop probabilities with lower probability than droptol relative
    % to the most likely state. 
    droptol = 1e-4;
    
    fprintf('Starting State-Space exploration\n');
    sse_timer = tic;
    last_percent_reported=0;
   
    p0 = speye(Ndofs,Ndofs);
    % Costruct matrix to be used in uniformization method
    lambda_max = max(abs(diag(D)));
    I = speye(Ndofs,Ndofs);
    A = I+D/lambda_max;
    nmax = 100;
    % tolp should be <= The DFSP tolerance \epsilon.  
    tolp=1e-2;
    pdv = p0;
    temp = p0;
    totp = 1;
    
    for i=0:nmax
        
        
       pp  = poisspdf(i,lambda_max*tau);
       temp = A*temp; 
       pdv = pdv+pp*temp;
       totp=totp-pp;
    
%        if(floor(i/nmax*100) > last_percent_reported + 4 )
%            %fprintf('%g%% complete\n',floor(i/Ndofs*100));
%            last_percent_reported = floor(i/nmax*100);
%            elapsed_time = toc(sse_timer);
%            time_remain = floor((1-(i/nmax))/((i/nmax)/elapsed_time));
%            fprintf('%g%% complete\t\t\telapsed: %g\teta: %is\n',last_percent_reported,elapsed_time,time_remain);
%        end    11
       
       if(totp<tolp)
           break;
       end
    end
    
    % Truncate 
    for i=1:Ndofs
        

       
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = urdme_init_DFSP_findM(fem,Ncells,Error_Tolerance,sumDiag,maxDiag,dt_in)
    % DFSP is faster than NSM if the expected number of diffusional 
    %  events (a_diff*tau_D) is greater than the number of voxels (K).
    %            a_diff * tau_D > K
    %            sumDiag * tau_D > K
    %            tau_D > K / sumDiag;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt =  max([ Ncells/sumDiag dt_in]);
    fprintf('urdme_init_DFSP_findM: dt = %g\n',dt);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    maximum_M = 25;
    minum_M = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find optimum MAX parameter
    M=-1;
    for i= minum_M:maximum_M
        tmp_err = 1-sum(poisspdf(0:i,maxDiag*dt));
        fprintf('\tM=%g E[err]=%g  (maxDiag=%g)\n',i,tmp_err,full(maxDiag));
        if(tmp_err < Error_Tolerance)
            M = i;
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
    dt_max  = min(diff(fem.urdme.tspan));
    %
    while(dt_set ~=1)
        if ~isfinite(dt)
            error('non-finite dt');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tmp_err = 1-sum(poisspdf(0:M,maxDiag*dt));
        fprintf('DFSP tau_D=%g\t\terror=%e\t\ttol=%e\n',dt,tmp_err,Error_Tolerance);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if( tmp_err < Error_Tolerance )
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
    fprintf('saving DFSP matrix as %s\n',dfsp_cache_filename);
    save(dfsp_cache_filename,'-STRUCT','fem');
end


