%  Assemble fvm discretization
%  
%
%  A first implementation. 
% 


function K = assemblefvm(umod,velo)

     
    [dim,ndofs] = size(umod.comsol.mesh.p);
       
    if dim == 2 
        K = assemblefvm2d(umod,velo);
    elseif dim == 3
        K = assemblefvm3d(umod,velo);
    else 
        error('assembly failed, the dimension must be 2 or 3.');
    end
        
            
        