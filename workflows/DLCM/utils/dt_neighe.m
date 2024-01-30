function [Ne,De] = dt_neighe(V,R)
%DT_NEIGHE Neighbor edge matrix.
%   Ne = DT_NEIGHE(V,R) returns a neighbor matrix the same sparsity
%   pattern as the neighbor matrix N, see DT_OPERATORS. For a given
%   population count vector u, v = Ne*u is the fraction of all voxel
%   boundaries occupied according to u. That is, v(i) is the populated
%   fraction of the total boundary of voxel i.
%
%   [Ne,De] = DT_NEIGHE(V,R) behaves differently near boundaries. In
%   this syntax, v = Ne*u+De has the same interpretation as above,
%   however, now all external boundaries are counted as if they were
%   populated. See the text output in parenthesis in the examples
%   below.
%
%   The inputs V and R are described in MESH2DUAL.
%
%   Example:
%     % circle geometry
%     C1 = [1 0 0 1]';
%     gd = [C1];
%     sf = 'C1';
%     ns = char('C1')';
%     G = decsg(gd,sf,ns);
%
%     % triangulation
%     [P,E,T] = initmesh(G,'hmax',0.5);
%     figure, pdegplot(G); hold on
%     triplot(T(1:3,:)',P(1,:),P(2,:),'k:');
%
%     % dual mesh
%     [V,R] = mesh2dual(P,E,T,'voronoi');
%     patch('Faces',R,'Vertices',V, ...
%           'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.25, ...
%           'EdgeColor','black');
%
%     % population vector
%     u = zeros(size(R,1),1); u(rand(size(u)) > 0.25) = 1;
%     patch('Faces',R(u == 1,:),'Vertices',V, ...
%           'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.25, ...
%           'EdgeColor','black');
%
%     % neighbor boundary cover
%     Ne = dt_neighe(V,R); [Ne_,De_] = dt_neighe(V,R);
%     v1 = Ne*u;
%     v2 = Ne_*u+De_;
%     for i = 1:numel(v1)
%       text(P(1,i),P(2,i),sprintf('%.2f',v1(i)));
%       text(P(1,i),P(2,i)-0.05,sprintf('(%.2f)',v2(i)));
%     end
%     axis tight, axis equal
%
%     % regular meshes
%     for mesh = 1:2
%       [P,E,T] = basic_mesh(mesh,5);
%       [V,R] = mesh2dual(P,E,T,'voronoi');
%       figure, patch('Faces',R,'Vertices',V, ...
%          'FaceColor',[0.9 0.9 0.9],'EdgeColor','black');
%       hold on,
%
%       u = zeros(size(R,1),1); u(rand(size(u)) > 0.25) = 1;
%       patch('Faces',R(u == 1,:),'Vertices',V, ...
%           'FaceColor',[0.2 0.2 0.2],'FaceAlpha',0.25, ...
%           'EdgeColor','black');
%       Ne = dt_neighe(V,R); [Ne_,De_] = dt_neighe(V,R);
%       v1 = Ne*u;
%       v2 = Ne_*u+De_;
%       for i = 1:numel(v1)
%         text(P(1,i),P(2,i),sprintf('%.2f',v1(i)));
%         text(P(1,i),P(2,i)-0.05,sprintf('(%.2f)',v2(i)));
%       end
%       axis([-1 1 -1 1]); axis equal
%     end

% S. Engblom 2023-11-13

% periodically wrap with NaN's
E = cat(3,R,[R(:,2:end) NaN(size(R,1),1)]);
% find all instances of a single NaN in 3rd dim
[i,j] = find(sum(isnan(E),3) == 1);
ix(i(end:-1:1)) = j(end:-1:1); % (look for smallest index)
% wrap around the first column into the position thus found
E((1:size(E,1))+size(E,1)*(ix-1)+prod(size(E,[1 2]))) = E(:,1,1);

% enforce a single orientation of all edges
E = sort(E,3);

% add information about patch #, permute & reshape and then drop all NaN's
E = cat(3,E,repmat((1:size(E,1))',1,size(E,2)));
E = reshape(permute(E,[3 2 1]),3,[]);
E(:,any(isnan(E))) = [];

% determine unique edges & compute their lengths
[U,ia,ib] = fsetop('unique',E(1:2,:));
edg = sqrt(sum(diff(reshape(V(U,:),2,[],2),[],1).^2,3));

% find all edges shared between patches
ix = find(E(3,:) ~= E(3,ia(ib)));
% [E(1:3,ix); E(3,ia(ib(ix)))] is a list of shared edges
% assemble sparse & symmetric edge neighbor matrix
Ne = sparse(E(3,ix),E(3,ia(ib(ix))),edg(ib(ix)),size(R,1),size(R,1));
Ne = Ne+Ne';

% output syntax
if nargout < 2
  % normalization excluding outer boundary edges:
  Ne = Ne./sum(Ne,2);
else
  % normalization by *total* edge length:
  De = full(diag(sparse(E(3,:),E(3,:),edg(ib))));
  Dn = full(sum(Ne,2));
  Ne = Ne./De;
  
  % vector of "constant" outer boundary cover:
  De = 1-Dn./De;
end
