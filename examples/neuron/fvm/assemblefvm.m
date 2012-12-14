%  Assemble fvm discretization
%  
%
%  A first implementation. 
% 


function K = assemblefvm(umod,velo)

    if iscmp4x(umod.comsol) 
      xmi = mphxmeshinfo(umod.comsol);
      [dim,ndofs] = size(xmi.dofs.coords);
    else
      [dim,ndofs] = size(umod.comsol.mesh.p);
    end
      
    if dim == 2 
        K = assemblefvm2d(umod,velo);
    elseif dim == 3
        K = assemblefvm3d(umod,velo);
    else 
        error('assembly failed, the dimension must be 2 or 3.');
    end
        
            
        