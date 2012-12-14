% URDME model file for the neuron example. 

function umod = neuron(umod)

global ft;
global cent;

Mspecies   = 3;
Mreactions = 8;

if iscmp4x(umod.comsol) %Comsol 4.x
  xmi = mphxmeshinfo(umod.comsol);
  dofs = xmi.dofs; 
  tri = xmi.elements.tri.nodes+1;
  nodes=dofs.nodes(1:Mspecies:end)+1;
else
  dofs = xmeshinfo(umod.comsol,'Out','dofs');
  tri = umod.comsol.mesh.e;
  nodes= dofs.nodes(1:Mspecies:end);
end
p = dofs.coords(:,1:Mspecies:end);

[dim,Ncells] = size(p); 
Ndofs = dim*Ncells;

V  = 1;
Vk = 2;
Vd = 3;

N = zeros(Mspecies,Mreactions);
N(V,1)=1;
N([V Vk],2) = [-1 1];
N([V Vd],3) = [-1 1];
N([V Vk],4) = [1 -1];
N([V Vd],5) = [1 -1];
N([Vk Vd],6)= [-1 1];
N([Vk Vd],7)= [1,-1];
N(V,8) = -1;

umod.N = sparse(N);
umod.G = sparse(ones(Mreactions,Mspecies+Mreactions));
    
% Subdomains
all=1:Ncells;

x = p(1,:);
y = p(2,:);
z = p(3,:);

umod.terminus = find(z<=0.7);
umod.axon     = find(z<1.15&(z>0.7&(x<-0.16&x>-0.3))&(z>0.7&(y>-2.1&y<-2.0)));
umod.axon     = union(umod.terminus,umod.axon);
umod.soma     = find((y>-2.17)&(y<-1.97)&(z>1.08)&(z<1.32)&(x>-0.3)&(x<-0.05));
dendrites=all;
dendrites = setdiff(dendrites,umod.axon);
umod.dendrites = setdiff(dendrites,umod.soma);
umod.soma  = setdiff(umod.soma,umod.axon); 
umod.data  = zeros(1,Ncells);
umod.data(umod.soma)=1;
umod.data(umod.axon)=2;

% Surface mesh data.
pn = 1:Ncells;
pn(nodes)=pn;
tri = pn(tri(1:3,:));
[~,ntri]=size(tri);

mem = find(umod.sd==2);

% To approximate a velocity field we construct tangent vectors to the
% surface (defined by the triangles) in directions determinied for 
% different regions of the neuron. Then for any other point
% we interpolate. 

% Boundary triangle representation
trep = TriRep(tri(1:3,:)',p');

%centroids 
cent = incenters(trep,(1:ntri)');

% Normals to surface triangle elements
fn = faceNormals(trep,(1:ntri)');

p0_soma = [-0.19,-2.05,1.21];
p0_axon = [10,-10,-10];

ft = zeros(size(fn));
for i=1:ntri
    
   if (any(umod.data(trep(i,:))==2))
     ft(i,:)=cross(cross(fn(i,:),p0_axon-cent(i,:)),fn(i,:)); 
   else
     ft(i,:)=cross(cross(fn(i,:),p0_soma-cent(i,:)),fn(i,:));
     % Smaller net speed in dendrites due to mixed polarity of fibers. 
     if(any(umod.data(trep(i,:))==0))
         ft(i,:)=0.01*ft(i,:);
     end
   end
   ft(i,:)=ft(i,:)/norm(ft(i,:));
   
end

% Delaunay triangulation object for nearest neighbour interpolation. 
global DT;
DT = DelaunayTri(cent);

% Active transport in axon 
addpath(['fvm']);

A1  = assemblefvm(umod,@vfk);
A2  = assemblefvm(umod,@vfd);

M  = spdiags(1./umod.vol,0,Ncells,Ncells); 
A1  = A1*M;
A2  = A2*M;

[i1,j1,s1]=find(A1);
i1 = Mspecies*(i1-1)+Vk;
j1 = Mspecies*(j1-1)+Vk;

[i2,j2,s2]=find(A2);
i2 = Mspecies*(i2-1)+Vd;
j2 = Mspecies*(j2-1)+Vd;

i = [i1;i2];
j = [j1;j2];
s = [s1;s2];

A = sparse(i,j,s,Ndofs,Ndofs);

umod.A = A;
umod.DM = umod.D;
umod.D = umod.D+umod.A;

% Initial condition
umod.u0 = zeros(3,Ncells);

% Microtubule density in axon.
mtd = 10;  
umod.data(2,umod.axon)=mtd;
umod.tspan = 0:1:100;

end

% Extremely simplified velocity field. 
function v = vfk(x)

global ft;
% Find direction vector by nearest neighbour interpolation in ft.
global DT;
P = nearestNeighbor(DT,x);
v = 0.1*[ft(P,1) ft(P,2) ft(P,3)]'; %% SIGN ???

end

% Extremely simplified velocity field. 
function v = vfd(x)

global ft;
% Find direction vector by nearest neighbour interpolation in ft.
global DT;
P = nearestNeighbor(DT,x);
v = -0.05*[ft(P,1) ft(P,2) ft(P,3)]'; %% SIGN ???

end
 


