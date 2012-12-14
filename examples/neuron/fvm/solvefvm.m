% FVM for simple convection problems (linear transport equation)


function fem = solvefvm(fem, TSPAN, velo)

%% Assembly
    
    if ~isfield(fem,'xmesh') || ~fem.xmesh.initialized
        fem.xmesh=meshextend(fem);
    end
    
    dofs  = xmeshinfo(fem,'Out','dofs');
    nodes = xmeshinfo(fem,'Out','nodes');
    Ndofs = length(dofs.nodes);    
    
    if ~isfield(fem,'urdme')
       fem.urdme = []; 
    end
    
    % Stiffness matrix
    if ~isfield(fem.urdme, 'K') 
        % Assemble (first order upwind)
         disp(['Assembling stiffness matrix...']);
         fem.urdme.K  = assemblefvm(fem,velo);   
         disp(['Done.']);
          
    end
    
    % Lumped mass matrix
    if ~isfield(fem.urdme, 'vol')
         disp(['Assembling mass matrix...']);
         % OBS! It is important that no stabilisation 
         % (such as streamline diffusion) is active in 
         % the FEM model. If it is, the mass (damping) matrix
         % will contain negative elements. 
         M = assemble(fem,'Out','D');
         fem.urdme.vol=sum(M,2);
         %fem.M  = spdiags(1./vol(nodes.dofs),0,Ndofs,Ndofs);
         disp(['done.']);
    end
   
    % Initial condition 
    u0 = asseminit(fem);
    % we are assembling in node ordering (as in fem.mesh.p), so we need to 
    % map u0 to this ordering. 
    u0 = u0.u(nodes.dofs);
  
    
%% Time integration. 

    disp(['Solving ...']);    
    
    M = spdiags(1./fem.urdme.vol(nodes.dofs),0,Ndofs,Ndofs);
    [t,U] = ode45(@fvmrhs,TSPAN,u0,[],M*fem.urdme.K);
    U = U';
    % Map back to COMSOL's dof-ordering for postprocessing. 
    fem.sol = femsol(U(dofs.nodes,:),'Tlist',t);
    
    