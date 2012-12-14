% ASSEMBLEFVM3D 
%
% Assemble the stiffness matrix for a first order upwind
% FVM discretization of the linear transport equation. 
% Neumann boundary conditions on all boundaries.
%

function A  = assemblefvm3d(fem,velo)



% Permutation to get the correct ordering of the dofs. 
if iscmp4x(fem.comsol)
  xmi = mphxmeshinfo(fem.comsol);
  dofs = xmi.dofs; 
  Mspecies = numel(xmi.dofs.dofnames);
  nodes=dofs.nodes(1:Mspecies:end)+1;
  t = xmi.elements.tet.nodes+1;
else
  dofs = xmeshinfo(fem.comsol,'Out','dofs');
  Mspecies = numel(dofs.names);
  nodes=dofs.nodes(1:Mspecies:end);
  t = umod.comsol.mesh.t;
end

p = dofs.coords(:,1:Mspecies:end);
[foo,Ncells]=size(p);
pn = 1:Ncells;
pn(nodes)=pn;
t = pn(t);

%Number of nodes, triangles and edges.
N = length(p);
T = length(t);

nt = size(t);

A   = sparse(N,N);
ML = zeros(4,4);
disp('Assembling active  transport...');
last_reported = 0;
for tri = 1:T

    % the vertices of tetrahedron tri
    tetra = t(1:4,tri);
    nodes = p(:,tetra);

    %the centroid of tetrahedron tri
    tetcenter=(1/4)*[sum(nodes(1,:));sum(nodes(2,:));sum(nodes(3,:))];

    % Assemble the 4x4 local stiffness matrix

    % Centroids of the facets of the primal tetrahedron.
    tcenter134 = (1/3)*[sum(nodes(1,[1 3 4]));sum(nodes(2,[1 3 4]));sum(nodes(3,[1 3 4]))];
    tcenter124 = (1/3)*[sum(nodes(1,[1 2 4]));sum(nodes(2,[1 2 4]));sum(nodes(3,[1 2 4]))];
    tcenter123 = (1/3)*[sum(nodes(1,[1 2 3]));sum(nodes(2,[1 2 3]));sum(nodes(3,[1 2 3]))];
    tcenter234 = (1/3)*[sum(nodes(1,[2 3 4]));sum(nodes(2,[2 3 4]));sum(nodes(3,[2 3 4]))];

    % Midpoint on the edges
    mid12 = 0.5*(nodes(:,1)+nodes(:,2));
    mid13 = 0.5*(nodes(:,1)+nodes(:,3));
    mid23 = 0.5*(nodes(:,2)+nodes(:,3));
    mid14 = 0.5*(nodes(:,1)+nodes(:,4));
    mid24 = 0.5*(nodes(:,2)+nodes(:,4));
    mid34 = 0.5*(nodes(:,3)+nodes(:,4));

    % The six facets of the dual voronoi element, ordered clockwise
    % with respect to the first vertex, i.e. facet12 with respect to
    % vertex 1 etc.
    facet12 = [mid12,tcenter123,tetcenter,tcenter124];
    facet13 = [mid13,tcenter134,tetcenter,tcenter123];
    facet23 = [mid23,tcenter123,tetcenter,tcenter234];
    facet14 = [mid14,tcenter124,tetcenter,tcenter134];
    facet24 = [mid24,tcenter234,tetcenter,tcenter124];
    facet34 = [mid34,tcenter134,tetcenter,tcenter234];

    % Areas and normals of facets of dual. With the ordering of nodes
    % in the facets as above, the normal is outward with respect to the
    % first vertex (as above).
    [n12,area12] = quadarea(facet12);
    [n13,area13] = quadarea(facet13);
    [n23,area23] = quadarea(facet23);
    [n14,area14] = quadarea(facet14);
    [n24,area24] = quadarea(facet24);
    [n34,area34] = quadarea(facet34);

    % Upwind discretization, velocity field "velo3d".
    % WHERE TO EVALUATE THE VELOCITY FIELD?
    
    v = velo([mid12 mid13 mid23 mid14 mid24 mid34]');
    
    flux = area12*n12*v(:,1);
    ML(2,1) = max(0,flux);
    ML(1,2) = max(0,-flux);
   
    flux = area13*n13*v(:,2);
    ML(3,1) = max(0,flux);
    ML(1,3) = max(0,-flux);
    
    flux = area23*n23*v(:,3);
    ML(3,2) = max(0,flux);
    ML(2,3) = max(0,-flux);
    
    flux = area14*n14*v(:,4);
    ML(4,1) = max(0,flux);
    ML(1,4) = max(0,-flux);
 
    flux = area24*n24*v(:,5);
    ML(4,2) = max(0,flux);
    ML(2,4) = max(0,-flux);

    flux = area34*n34*v(:,6);
    ML(4,3) = max(0,flux);
    ML(3,4) = max(0,-flux);
  
    A(tetra,tetra) = A(tetra,tetra)+ML;
    
    completed = floor(100*tri/T);
    if completed > last_reported+1;
     disp([num2str(completed) ' % completed.']);
     last_reported = completed;
    end
     
end

A = A-spdiags(sum(A)',0,N,N);
disp('completed.')
