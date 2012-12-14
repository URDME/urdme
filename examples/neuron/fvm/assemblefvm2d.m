%  Assemble fvm discretization
%  
%
%  A first implementation. 
% 


function A = assemblefvm2d(fem,velo)


% Permutation to get the correct ordering of the dofs. 
dofs = xmeshinfo(fem,'out','dofs');
Mspecies = numel(dofs.names);
p = dofs.coords(:,1:Mspecies:end);
[~,Ncells]=size(p);
pn = 1:Ncells;
pn(dofs.nodes(1:Mspecies:end))=pn;
t = pn(fem.mesh.t);

%Number of nodes, triangles and edges.
N = length(p);
T = length(t);

%A = zeros(N,N);
nt = size(t);

% No spalloc in COMSOL ??
A   = spalloc(N,N,8*N);

% Local node-matrix
ML = zeros(3,3);
   
%bar = waitbar(0,'Assembling active  transport...');
for tri = 1:T
   
   % the vertices of triangle tri
   triangle = t(1:3,tri);
   nodes = p(:,triangle);
   
   %the centroid of triangle tri
   tcenter=(1/3)*[sum(nodes(1,:));sum(nodes(2,:))];
   
   % Assemble the 3x3 local stiffness matrix
   ML(:,:) = 0;
   % Midpoint on edge ij.
   mid12 = 0.5*(nodes(:,1)+nodes(:,2));
   mid13 = 0.5*(nodes(:,1)+nodes(:,3));
   mid23 = 0.5*(nodes(:,2)+nodes(:,3));
   
   % edge vectors (tangents)
   edge12 = mid12-tcenter;
   edge13 = mid13-tcenter;
   edge23 = mid23-tcenter;
   
   % outwards normals
   n12 = [-edge12(2),edge12(1)];
   n12=n12/norm(n12);
   if dot(n12,mid12-nodes(:,1)) < 0
       n12=-n12;
   end
   
   n13 = [-edge13(2),edge13(1)]; 
   n13=n13/norm(n13);
   if dot(n13,mid13-nodes(:,1)) < 0
       n13=-n13;
   end
   
   n23 = [-edge23(2),edge23(1)]; 
   n23=n23/norm(n23);
   if dot(n23,mid23-nodes(:,2)) < 0
       n23=-n23;
   end
   
   % distance between centroid and midpoints    
   dmc12 = norm(edge12); 
   dmc13 = norm(edge13); 
   dmc23 = norm(edge23); 
   
   % Upwind discretization, velocity field "velo".
   flux = n12*velo(mid12);
   if (flux >= 0)
       ML(2,1)=dmc12*flux;
   else
       ML(1,2)=-dmc12*flux;
   end

   flux = n13*velo(mid13);
   if (flux >= 0)
       ML(3,1)=dmc13*flux;
   else
       ML(1,3)=-dmc13*flux;
   end

   flux = n23*velo(mid23);
   if (flux >= 0)
       ML(3,2)=dmc23*flux;
   else
       ML(2,3)=-dmc23*flux;
   end
  
   % Assemble into global matrix
   A(triangle,triangle) = A(triangle,triangle)+ML;
   %waitbar(tri/T,bar);
 
end

A = A-spdiags(sum(A)',0,N,N);




    






