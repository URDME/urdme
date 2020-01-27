function [V,R] = mesh2dual(P,E,T,type)
%MESH2DUAL Dual mesh from triangulation.
%   [V,R] = MESH2DUAL(P,E,T,Type) constructs the dual mesh (V,R) from
%   the PDE Toolbox format triangulation (P,E,T). Use E = zeros(7,0)
%   if the edge matrix is missing, in which case boundary triangles
%   will not be handled.
%
%   The dual mesh follows the PATCH-function format: V is an NV-by-2
%   vertex matrix containing the x- and y-coordinates as columns. R is
%   an NR-by-MaxElem matrix of mesh elements, where each row R(i,:)
%   points into V thus forming element i. The rows are padded with
%   NaNs as needed.
%
%   The input Type is either 'median' (the default) or 'voronoi'. Type
%   == 'median' constructs the dual by joining the midpoints of the
%   edges with triangle median points, adding also nodes at the
%   boundary whenever required. This type of mesh is standard in
%   Finite Volume Methods.
%
%   Type == 'voronoi' joins instead the midpoints of the circles
%   circumscribing the triangles, but taking extra care at the edges
%   of the triangulation. At edges, whenever a midpoint is on the
%   outside of the boundary, the midpoint is projected back
%   orthogonally towards the two triangle edges on the inside of the
%   domain, and the two edge intersection points thus obtained are
%   used instead. Note: it is not guaranteed that no dual mesh element
%   stretches outside the boundary since non-boundary triangles of
%   very low quality can produce such elements; only boundary
%   triangles are handled carefully.
%
%   For both types of dual meshes the ordering in the mesh is
%   preserved in the sense that the midpoint of element i is P(:,i).
%
%   Examples:
%     % geometry
%     C1 = [1 0 0 1]';
%     C2 = [1 -0.4 0 0.25]';
%     C3 = [1 0.4 0 0.25]';
%     gd = [C1 C2 C3];
%     sf = '(C1+C2)-C3';
%     ns = char('C1','C2','C3')';
%     G = decsg(gd,sf,ns);
%
%     % triangulation
%     [P,E,T] = initmesh(G,'hmax',0.15);
%     figure, clf, pdegplot(G); hold on
%     triplot(T(1:3,:)',P(1,:),P(2,:),'k:');
%
%     % dual mesh
%     [V,R] = mesh2dual(P,E,T,'voronoi');
%     patch('Faces',R,'Vertices',V, ...
%           'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.25, ...
%           'EdgeColor','black');
%     axis tight, axis equal
%
%   See also INITMESH, PATCH, VORONOI, DECSG.

% S. Engblom 2017-10-22 (Very major revision)
% S. Engblom 2008-10-07

% mesh tolerance
TOL = eps(1024)*max(abs(P),[],2);

% subdomain info is ignored:
T = T(1:3,:);

% type median/voronoi
if nargin < 4
  type = 1;
else
  switch lower(type)
   case 'median', type = 1;
   case 'voronoi', type = 2;
   otherwise, error('Unknown dual mesh type.');
  end
end

% Voronoi: precompute midpoints of all circumscribed circles
if type == 2
  AB = P(:,T(2,:))-P(:,T(1,:));
  nB = sum(AB.^2,1);
  AC = P(:,T(3,:))-P(:,T(1,:));
  nC = sum(AC.^2,1);
  D = 2*(AB(1,:).*AC(2,:)-AB(2,:).*AC(1,:));
  U = [AC(2,:).*nB-AB(2,:).*nC; ...
       AB(1,:).*nC-AC(1,:).*nB];
  U = tprod(U,1./D,[1 2],[3 2]);
  Tmid = U+P(:,T(1,:));
  % formulas taken from
  % https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_2

  % the rest of the code now deals with midpoints outside boundary
  % triangles...

  % find all triangles with the midpoint outside the triangle (i.e.,
  % this is simply a condition of obtuse triangles)
  d = sum((P(:,T(2,:))-P(:,T(3,:))).^2,1);
  [d,id] = sort([nB; nC; d]);
  iout = find(sum(d(1:2,:),1) < d(3,:));  % (Pythogarean inequality)

  % determine the edges where the midpoints were outside (longest edges)
  edg = [2 1; 3 1; 2 3]'; % (ordering as d above was formed)
  long = sort(T(tsum(edg(:,id(3,iout)),(iout-1)*3,[1 2],[3 2])));

  % intersect those edges with all boundary edges...
  [long,it,ie] = fsetop('intersect',long,sort(E(1:2,:)));
  iout = iout(it);
  % ...these are the problematic ones dealt with now:

  % prepare a (possibly empty) lookup table (Tout,PT)
  Tout = sparse(1,iout,1:numel(iout),1,size(T,2));
  PT = zeros(5,0);
  % Tout is 1-by-NT, where Tout(i) is ix ~= 0 if the ith triangle is
  % outside a boundary edge, in which case PT(1,ix) contains the edge
  % number, and P(2:3,ix) and P(4:5,ix) the two projected edge
  % intersection points
  if ~isempty(iout)
    % for the list of problematic triangles, additionally compute the two
    % points of intersection with the longest edge starting from the
    % midpoint and towards the midpoints of the two shorter edges
    short = T(tsum(reshape(edg(:,id(1:2,iout)),2,2,[]), ...
                   (iout-1)*3,[1 2 3],[4 3]));
    % midpoints of those short edges:
    short = reshape(mean(reshape(P(:,short),2,2,[]),2),2,2,[]);

    % now solve Tmid+alpha1*v1 = P1+beta1*w with w = P2-P1 the longe edge
    % and (v1,v2) the midpoints of the two short edges
    v1 = reshape(short(:,1,:),2,[])-Tmid(:,iout);
    v2 = reshape(short(:,2,:),2,[])-Tmid(:,iout);
    P1 = P(:,long(1,:));
    w = P(:,long(2,:))-P1;

    % 2-by-2 inversion
    W = permute(cat(3,w,-v1),[1 3 2]);
    invW = tprod(1./(W(1,1,:).*W(2,2,:)-W(2,1,:).*W(1,2,:)), ...
                 [W(2,2,:) -W(1,2,:); -W(2,1,:) W(1,1,:)], ...
                 [4 5 3],[1 2 3]);
    ab = tprod(invW,Tmid(:,iout)-P1,[1 -1 2],[-1 2]);
    beta1 = ab(1,:);

    W = permute(cat(3,w,-v2),[1 3 2]);
    invW = tprod(1./(W(1,1,:).*W(2,2,:)-W(2,1,:).*W(1,2,:)), ...
                 [W(2,2,:) -W(1,2,:); -W(2,1,:) W(1,1,:)], ...
                 [4 5 3],[1 2 3]);
    ab = tprod(invW,Tmid(:,iout)-P1,[1 -1 2],[-1 2]);
    beta2 = ab(1,:);

    % crossing points thus computed go into the table
    PT = [ie; ...
          P1+tprod(beta1,w,[3 2],[1 2]); ...
          P1+tprod(beta2,w,[3 2],[1 2])];
  end
end

% list of dual elements
ME = cell(1,size(P,2));

% one dual element per node in the primal mesh, but the nodes on the
% edges towards subdomain 0 need to be handled differently
ip_ = 1:size(P,2);

% loop over all nodes on edges to subdomain 0
for ie = find(any(E(6:7,:) == 0,1))
  ip = E(1,ie);
  it = find(any(T == ip,1));
  if type == 1
    % unique nodes in the union of all triangles sharing this node
    j = fsetop('unique',reshape(T(:,it),1,[]));
    j = j(j ~= ip);

    % to sort correctly, add an extra point outside the edge
    v = P(:,E(2,ie))-P(:,ip); % v along this edge
    w = [v(2); -v(1)]; % w orthogonal to this edge
    if E(7,ie) ~= 0
      w = -w; % ensure that w points towards subdomain 0
    end
    % one tiny fraction _outwards_ of edge
    z = 0.5*v+TOL.*w/norm(w);

    % put this point first, then the rest of the nodes
    z = [complex(z(1),z(2)) complex(P(1,j)-P(1,ip),P(2,j)-P(2,ip))];
    % sort nodes by angle
    [~,is] = sort(angle(z));
    % wrap where the extra point went
    i1 = find(is == 1);
    is = is([i1+1:end 1:i1-1])-1;
    % this is now a contiguous ordering of all edges:
    j = j(is);

    % edge midpoints, triangle centroids, and node itself
    za = 1/2*tsum(P(:,j),P(:,ip),[1 2],[1 3]);
    zb = 1/3*tsum(P(:,j(1:end-1))+P(:,j(2:end)),P(:,ip),[1 2],[1 3]);
    zz = reshape([za; [zb P(:,ip)]],2,[]);
  else
    % Voronoi: midpoints of circumscribed circles
    zz = Tmid(:,it);

    % are any of those midpoints outside of triangle and domain?
    [~,ii,ix] = find(Tout(it));
    if ~isempty(ii)
      % if so, use the points in the table and select the point closest to
      % the current node instead
      PT_ = reshape(PT(2:end,ix),2,2,[]);
      dist = sum(tsum(PT_,-P(:,ip),[1 2 3],[1 4]).^2,1);
      [~,jj] = min(dist,[],2);
      zz(:,ii) = PT_(:,jj(:)'+2*(0:numel(ix)-1));
    end
    je = PT(1,ix); % prevent unnecessary nodes

    % add node itself and (possibly) midpoints of boundary edges
    zz = [P(:,ip) zz];
    if all(ie ~= je)
      inext = E(2,ie);  
      zz = [zz 0.5*(P(:,ip)+P(:,inext))];
    end
    ibefore = find(E(2,:) == ip);
    if all(ibefore ~= je)
      ibefore = E(1,ibefore);
      zz = [zz 0.5*(P(:,ip)+P(:,ibefore))];
    end

    % sort nodes by their angle around their common mean
    [~,is] = sort(atan2(zz(2,:)-mean(zz(2,:)),zz(1,:)-mean(zz(1,:))));
    zz = zz(:,is);
  end

  % form macroelement
  ME{ip} = zz;

  % remove this node from the mesh
  ip_(ip) = 0;
end

% now loop over the remaining nodes in the primal mesh
ip_ = ip_(ip_ ~= 0);
for ip = ip_
  it = find(any(T == ip,1));
  if type == 1
    % unique nodes in the union of all triangles sharing this node
    j = fsetop('unique',reshape(T(:,it),1,[]));
    j = j(j ~= ip);
    % sort nodes by angle
    [~,is] = sort(atan2(P(2,j)-P(2,ip),P(1,j)-P(1,ip)));
    j = j(is);

    % edge midpoints and triangle centroids
    za = 1/2*tsum(P(:,j),P(:,ip),[1 2],[1 3]);
    zb = 1/3*tsum(P(:,j)+P(:,j([2:end 1])),P(:,ip),[1 2],[1 3]);
    zz = reshape([za; zb],2,[]);
  else
    % Voronoi: midpoints of circumscribed circles
    zz = Tmid(:,it);

    % are any of those midpoints outside of triangle and domain?
    [~,ii,ix] = find(Tout(it));
    if ~isempty(ii)
      % if so, use the two points in the table instead
      zz(:,ii) = [];
      zz = [zz reshape(PT(2:end,ix,:),2,[])];
    end
    % sort nodes by angle
    [~,is] = sort(atan2(zz(2,:)-P(2,ip),zz(1,:)-P(1,ip)));
    zz = zz(:,is);
  end

  % form macroelement
  ME{ip} = zz;
end

% produce (V,R)-format
elem_siz = cellfun('size',ME,2);

% unique vertices
ME = cat(2,ME{:});
[~,V_,R_] = fsetop('unique',fix(tprod(ME,1./TOL,1:2,1)));
V = ME(:,V_)';

% follow the patch-function's "compact format"
R = nan(size(P,2),max(elem_siz));
ir = 1;
for i = 1:size(R,1)
  rsiz = elem_siz(i);
  R(i,1:rsiz) = R_(ir:ir+rsiz-1);
  ir = ir+rsiz;
end

return;

% midpoint of circumscribed circle formulas at work
A = [0 0];
B = [1 0];
C = [0.5 0.15];

BB = B-A;
CC = C-A;
D = 2*(BB(1)*CC(2)-BB(2)*CC(1));
U = [CC(2)*norm(BB)^2-BB(2)*norm(CC)^2 ...
     BB(1)*norm(CC)^2-CC(1)*norm(BB)^2]/D;
zz = U+A;
r = norm(U);

theta = linspace(0,2*pi);
p = complex(zz(1),zz(2))+r*exp(1i*theta);
figure, plot(A(1),A(2),'bo'); hold on,
plot(B(1),B(2),'bo'); plot(C(1),C(2),'bo');
plot(zz(1),zz(2),'ro');
plot(p,'b--');
axis equal

% an example for which the Voronoi mesh is somewhat surprising:
% 6 obtuse triangles (2 serious ones in the middle):
P = [0 1; 0.1 0; 0 -1; -0.1 0; ...
     -1 0; 1 0; ...
     -1 1; 1 1; -1 -1; 1 -1; ...
     0 -2; 0 2]';
E = [11 9; 9 5; 5 7; 7 12; 12 8; 8 6; 6 10; 10 11]';
E = [E; zeros(1,size(E,2)); ones(2,size(E,2)); ...
     zeros(1,size(E,2)); ones(1,size(E,2))];
T = [1 2 3; 1 3 4; 1 2 6; 2 3 6; 1 5 4; 4 3 5; ...
     5 7 1; 1 6 8; 6 10 3; 3 9 5; ...
     3 11 9; 3 10 11; 12 8 1; 1 7 12]';
figure, clf,
triplot(T(1:3,:)',P(1,:),P(2,:),'k:'); hold on
[V,R] = mesh2dual(P,E,T,'voronoi');
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.25, ...
      'EdgeColor','black');
axis([-2 2 -3 3]);
