function flag = contact(P0,R0,PR0,Theta0,P1,R1,PR1,Theta1,dTheta,debug)
%CONTACT Determine how two cells interact.
%   FLAG = CONTACT(P0,R0,PR0,Theta0,P1,R1,PR1,Theta1,dTheta) returns
%     0: if cell 1 does not interact with cell 0,
%     1: if cell 1 interacts with cell 0 via protrusions,
%     2: if cell 1 interacts with cell 0 via junctional contact,
%     3: if cell 1 interacts with cell 0 via both protrusions and
%     junctional contacs.
%
%   A protrusional contact (FLAG = 1 or 3) means that either the
%   protrusions touch symmetrically or the protrusions of cell 1
%   touches the membrane of cell 0.
%
%   A junctional contact (FLAG = 2 or 3) means that either the
%   membrane of the cells touch symmetrically or the membrane of cell
%   1 touches the protrusions of cell 0.
%
%   The first cell is centered at point P0, with radius R0 and with a
%   protrusion of radius PR0 in the direction Theta0, understood also
%   to protrude in the inverse direction Theta0-pi; similarly for the
%   second cell (P1,R1,PR1,Theta1). The angular width of the
%   protrusions is dTheta.
%
%   This is a utility function and no error-checking is performed.

% S. Engblom 2016-11-22

% quick return
d = norm(P0-P1);
% no protrusion can possibly reach the other
if d > PR0+PR1
  flag = 0;
  return;
end

% hard-wired case
if dTheta >= pi
  if d <= PR0+R1
    flag = 3; % (= 1+2 since this contact will also be protrusional)
  elseif  d <= R0+PR1 || d <= PR0+PR1
    flag = 1;
  end
  return;
end

% assume no contact
flag = 0;

% discretize the angle dTheta using a minimum angle minangle
minangle = pi/12;
nangles = ceil(dTheta/minangle)+1;

% (1) We now solve the equation P0+C0*W0 = P1+C1*W1 for W0 and W1 all
% protrusional rays according to the above discretization, and hence
% find nangles^2 different pairs (C0,C1). This implies the same number
% of protrusional ray crossing points and if a single one of them is
% within the protrusion radius of both cells simultaneously, a
% protrusional contact is made.

% start by building a matrix with columns orthogonal to W0
theta0 = linspace(Theta0,Theta0+dTheta,nangles);
orthW0 = frepmat([sin(theta0); -cos(theta0)],[1 nangles]);

% next the W1 rays
theta1 = linspace(Theta1,Theta1+dTheta,nangles);
W1 = reshape(frepmat([cos(theta1); sin(theta1)],nangles),2,[]);

% solve by using the orthogonality of orthW0
C1 = (orthW0'*(P0-P1))./tprod(orthW0,W1,[-1 1],[-1 1]);

% the nangles^2 crossing points
pts = tsum(P1,tprod(W1,C1,[1 2],[2]),[1],[1 2]);
% the 2*nangles^2 distances to centers
dist = reshape(sqrt(sum(tsum(pts,-[P0 P1],[1 2],[1 3]).^2)),[],2);

% criteria
if any(dist(:,1) <= PR0 & dist(:,2) <= PR1)
  % protrusion touches symmetrically
  flag = 1;
else
  % (2) Otherwise we need to investigate if the protrusion of 1 touches
  % the membrane of 0. Determine the closest point to cell 0 in the
  % protrusion region of cell 1; if this point is within the radius of
  % cell 0, the cells are in protrusional contact.

  % solve min_c1 norm(P1+c1*W1-P0), or -c1 = Proj_W1 (P1-P0)
  W1 = W1(:,1:nangles:end);
  c1 = -W1'*(P1-P0);
  ibnd = find(abs(c1) > PR1);
  c1(ibnd) = PR1*sign(c1(ibnd));
  % the nangles closest points
  pts1 = tsum(P1,tprod(W1,c1,[1 2],[2]),[1],[1 2]);
  dist1 = reshape(sqrt(sum(tsum(pts1,-P0,[1 2],[1 3]).^2)),[],1);

  % criteria
  if any(dist1 <= R0)
    flag = 1;
  end
end

% next investigate junctional contacts
if d <= R0+R1
  % immediate case, membrane touches symmetrically
  flag = flag+2;
else
  % (3) We now need to investigate if the protrusion of 0 touches the
  % membrane of 1. This is done as in (2) above.

  % solve min_c0 norm(P0+c0*W0-P1), or -c0 = Proj_W0 (P0-P1)
  W0 = [-orthW0(2,:); orthW0(1,:)];
  W0 = W0(:,1:nangles);
  c0 = W0'*(P1-P0);
  ibnd = find(abs(c0) > PR0);
  c0(ibnd) = PR0*sign(c0(ibnd));
  pts0 = tsum(P0,tprod(W0,c0,[1 2],[2]),[1],[1 2]);
  dist0 = reshape(sqrt(sum(tsum(pts0,-P1,[1 2],[1 3]).^2)),[],1);

  % criteria
  if any(dist0 <= R1)
    flag = flag+2;
  end
end

if nargin < 10 || ~debug
  return;
end

% visualization under debug
w = [cos(Theta0) cos(Theta0+dTheta) cos(-pi+Theta0) cos(-pi+Theta0+dTheta); ...
     sin(Theta0) sin(Theta0+dTheta) sin(-pi+Theta0) sin(-pi+Theta0+dTheta)];
hold on
uv = PR0*w(:,1:2);
quiver(P0([1 1])',P0([2 2])',uv(1,:),uv(2,:),0,'b','LineWidth',2);
if debug == 2
  plot(P0(1)+PR0*cos(Theta0+linspace(0,dTheta,4*nangles)), ...
       P0(2)+PR0*sin(Theta0+linspace(0,dTheta,4*nangles)), ...
       'k','LineWidth',1);
  plot([P0(1) P0(1)+PR0],[P0(2) P0(2)],'k--','LineWidth',1);
  plot(P0(1)+0.5*PR0*cos(linspace(0,Theta0,20)), ...
       P0(2)+0.5*PR0*sin(linspace(0,Theta0,20)), ...
       'k','LineWidth',1);
end
uv = PR0*w(:,3:4);
quiver(P0([1 1])',P0([2 2])',uv(1,:),uv(2,:),0,'b','LineWidth',2);
% $$$ plot(P0(1)+PR0*cos(-pi+Theta0+linspace(0,dTheta,4*nangles)), ...
% $$$      P0(2)+PR0*sin(-pi+Theta0+linspace(0,dTheta,4*nangles)), ...
% $$$      'k','LineWidth',2);

w = [cos(Theta1) cos(Theta1+dTheta) cos(-pi+Theta1) cos(-pi+Theta1+dTheta); ...
     sin(Theta1) sin(Theta1+dTheta) sin(-pi+Theta1) sin(-pi+Theta1+dTheta)];
uv = PR1*w(:,1:2);
quiver(P1([1 1])',P1([2 2])',uv(1,:),uv(2,:),0,'r','LineWidth',2);
% $$$ plot(P1(1)+PR1*cos(Theta1+linspace(0,dTheta,4*nangles)), ...
% $$$      P1(2)+PR1*sin(Theta1+linspace(0,dTheta,4*nangles)), ...
% $$$      'k','LineWidth',2);
uv = PR1*w(:,3:4);
quiver(P1([1 1])',P1([2 2])',uv(1,:),uv(2,:),0,'r','LineWidth',2);
% $$$ plot(P1(1)+PR1*cos(-pi+Theta1+linspace(0,dTheta,4*nangles)), ...
% $$$      P1(2)+PR1*sin(-pi+Theta1+linspace(0,dTheta,4*nangles)), ...
% $$$      'k','LineWidth',2);

return;

% protrusion crossings
ii = find(R0 <= dist(:,1) & dist(:,1) <= PR0 & ...
          R1 <= dist(:,2) & dist(:,2) <= PR1);
plot(pts(1,ii),pts(2,ii),'ko','LineWidth',2,'MarkerSize',6)

if exist('pts1','var')
  ii = find(dist1 <= R0);
  plot(pts1(1,ii),pts1(2,ii),'rx','LineWidth',2,'MarkerSize',6);
end

if exist('pts0','var')
  ii = find(dist0 <= R1);
  plot(pts0(1,ii),pts0(2,ii),'b+','LineWidth',2,'MarkerSize',6);
end
