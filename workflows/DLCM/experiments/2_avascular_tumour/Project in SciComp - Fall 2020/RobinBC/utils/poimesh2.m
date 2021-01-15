function [p,e,t]=poimesh2(g,n1,n2)
%POIMESH Make regular mesh on a rectangular geometry.
%
%       [P,E,T]=POIMESH(G,NX,NY) constructs a regular mesh
%       on the rectangular geometry specified by G, by dividing
%       the "x edge" into NX pieces and the "y edge" into NY pieces,
%       and placing (NX+1)*(NY+1) points at the intersections.
%
%       The x edge is the one that makes the smallest angle
%       with the x axis.
%
%       [P,E,T]=POIMESH(G,N) uses NX=NY=N, and
%       [P,E,T]=POIMESH(G) uses NX=NY=1;
%
%       For best performance with POISOLV, the larger of NX and NY
%       should be a power of 2.
%
%       If G does not seem to describe a rectangle, P is zero on return.
%
%       See also INITMESH, POISOLV

%       Copyright 1994-2016 The MathWorks, Inc.

if nargin==2
  n2=n1;
elseif nargin==1
  n1=1;
  n2=1;
end

nbs=pdeigeom(g); % Number of boundary segments
% A rectangle has 4 BS
if nbs~=4
  p=0;
  return
end

d=pdeigeom(g,1:nbs);
% SD 0 must be on one side
if min(d(3:4,:))~=zeros(1,nbs)
  p=0;
  return
end
% and SD 1 on the other
if max(d(3:4,:))~=ones(1,nbs)
  p=0;
  return
end

% Start to identify points
[x,y]=pdeigeom(g,[1:nbs;1:nbs],d(1:2,:));
small=1000*eps; % Equality tolerance
scale=max(max(max(abs(x))),max(max(abs(y))));
tol=small*scale;

p=[x(1,1);y(1,1)]; % A first point
np=1;

for i=1:nbs
  for j=1:2
    xy=[x(j,i);y(j,i)];
    l=find(sum(abs(p(1:2,:)-xy*ones(1,size(p,2))))<tol);
    if isempty(l)
      p=[p [xy]]; % New point
      np=np+1;
      d(4+j,i)=np;
    else
      d(4+j,i)=l; % Old point
    end
  end
end
% Row 5 and 6 of d now contains indices of starting and ending points
% of each BS

if np~=4
  p=0;
  return
end

% Make SD 1 be on the left side
ix=find(d(4,:));
d([2 1 4 3 6 5],ix)=d(1:6,ix);

ev=[p(:,d(6,:))-p(:,d(5,:))];

% Small check for line-ness
nn=10; % for example
[x,y]=pdeigeom(g,(1:4)'*ones(1,nn), ...
     d(2,:)'*linspace(0,1,nn)+d(1,:)'*linspace(1,0,nn));
for k=1:4
  ev1=[diff(x(k,:));diff(y(k,:))];
  if any(ev1-ev(:,k)*ones(1,nn-1)/(nn-1)>small)
    p=0;
    return
  end
end

ev1=ev./(ones(2,1)*sqrt(ev(1,:).^2+ev(2,:).^2));

i1=find(ev1(1,:)==max(ev1(1,:)));
i1=i1(1);

% Find cycle
i2=find(d(5,:)==d(6,i1));
i3=find(d(5,:)==d(6,i2));
i4=find(d(5,:)==d(6,i3));

% last check for rectangle-ness
if abs(ev1(:,i2)'*ev1(:,i1))>small
  p=0;
  return
end
if abs(ev1(:,i3)'*ev1(:,i1)+1)>small
  p=0;
  return
end
if abs(ev1(:,i4)'*ev1(:,i1))>small
  p=0;
  return
end

p1=p(:,d(5,[i1 i2 i3 i4]));

% Generate mesh
np=(n1+1)*(n2+1);
ne=2*(n1+n2);
nt=2*n1*n2;

p=zeros(2,np);
e=zeros(7,ne);
t=zeros(4,nt);

ip=0:(np-1);
ip2=floor(ip/(n1+1));
ip1=(ip-(n1+1)*ip2);
p=p1(:,1)*ones(1,np)+[p1(:,2)-p1(:,1) p1(:,3)-p1(:,2)]*[ip1/n1; ip2/n2];

e(:,1:n1)= ...
[1:n1;2:(n1+1); ...
d(1,i1)+(0:(n1-1))/n1*(d(2,i1)-d(1,i1)); ...
d(1,i1)+(1:n1)/n1*(d(2,i1)-d(1,i1)); ...
i1*ones(1,n1);ones(1,n1);zeros(1,n1)];
e(:,n1+1:n1+n2)= ...
[n1+1:n1+1:n2*(n1+1);2*(n1+1):n1+1:(n2+1)*(n1+1); ...
d(1,i2)+(0:(n2-1))/n2*(d(2,i2)-d(1,i2)); ...
d(1,i2)+(1:n2)/n2*(d(2,i2)-d(1,i2)); ...
i2*ones(1,n2);ones(1,n2);zeros(1,n2)];
e(:,n1+n2+1:2*n1+n2)= ...
[(n2+1)*(n1+1):-1:n2*(n1+1)+2;(n2+1)*(n1+1)-1:-1:n2*(n1+1)+1; ...
d(1,i3)+(0:(n1-1))/n1*(d(2,i3)-d(1,i3)); ...
d(1,i3)+(1:n1)/n1*(d(2,i3)-d(1,i3)); ...
i3*ones(1,n1);ones(1,n1);zeros(1,n1)];
e(:,2*n1+n2+1:ne)= ...
[n2*(n1+1)+1:-(n1+1):n1+2;(n2-1)*(n1+1)+1:-(n1+1):1; ...
d(1,i4)+(0:(n2-1))/n2*(d(2,i4)-d(1,i4)); ...
d(1,i4)+(1:n2)/n2*(d(2,i4)-d(1,i4)); ...
i4*ones(1,n2);ones(1,n2);zeros(1,n2)];

% Change the connectivity matrix such that triangulation goes in the 
% other direction
it=0:(nt/2-1);
it2=floor(it/n1);
it1=(it-n1*it2);

% Mesh generation depends on Nvoxels being odd/even

if mod(n1+1,2) % Nvoxels is odd

    % Divide triangles into intervals
    start_int = (1:n1:nt); 
    
    % Start indices of intervals that starts with an odd number to keep
    odd_start = start_int(1:2:end);
    add_on_odd = 0:2:(n1-1); % Keep every other odd number

    % Define which odd indices to keep in each interval
    matr1 = odd_start'*ones(1, length(add_on_odd));
    keep_odd = matr1 + add_on_odd;
    keep_odd = keep_odd(:);

    % Start indices of intervals that starts with an even number to keep
    even_start = start_int(2:2:end);
    add_on_even = 1:2:n1;
    
    % Define which even indices to keep in each interval
    matr2 = even_start'*ones(1, length(add_on_even));
    keep_even = matr2 + add_on_even;
    keep_even= keep_even(:);

    % Final definition of which triangles to change 
    keep = union(keep_odd, keep_even);
    change = setdiff((1:nt), keep);

    % Divide into different columns for upper and lower triangles
    keep_low = keep(keep <= nt/2);
    keep_up = keep(keep > nt/2);

    change_low = change(change <= nt/2);
    change_up = change(change > nt/2);

    % Lower triangles /| and |\
    t(1:4,keep_low) = [it1(keep_low) + (n1+1)*it2(keep_low) + 1; ...
                       it1(keep_low) + (n1+1)*it2(keep_low) + 2; ...
                       it1(keep_low) + (n1+1)*(it2(keep_low)+1) + 2; ...
                       ones(1,length(keep_low))];

    t(1:4,change_low) = [it1(change_low) + (n1+1)*it2(change_low) + 1; ...
                         it1(change_low) + (n1+1)*it2(change_low) + 2; ...
                         it1(change_low) + (n1+1)*(it2(change_low)+1) + 1; ...
                         ones(1,length(change_low))];

    % Upper triangles |/ and \|
    t(1:4,keep_up) = [it1(keep_low) + (n1+1)*it2(keep_low) + 1; ... 
                      it1(keep_low) + (n1+1)*(it2(keep_low)+1) + 2; ...
                      it1(keep_low) + (n1+1)*(it2(keep_low)+1) + 1; ...
                      ones(1,length(keep_low))];

    t(1:4,change_up) = [it1(change_low) + (n1+1)*it2(change_low) + 2; ... 
                        it1(change_low) + (n1+1)*(it2(change_low)+1) + 2; ...
                        it1(change_low) + (n1+1)*(it2(change_low)+1) + 1; ...
                        ones(1,length(change_low))];


else % Nvoxels is even
    end_column1 = (nt/2+2):2:nt;
    end_column2 = 2:2:(nt/2);
    
    % Lower triangles /| and |\
    t(1:4,1:2:(nt/2))=[it1(1:2:(nt/2)) + (n1+1)*it2(1:2:(nt/2)) + 1; ...
                       it1(1:2:(nt/2)) + (n1+1)*it2(1:2:(nt/2)) + 2; ...
                       it1(1:2:(nt/2)) + (n1+1)*(it2(1:2:(nt/2))+1) + 2; ...
                       ones(1,length(1:2:(nt/2)))];

    t(1:4,end_column2)=[it1(end_column2) + (n1+1)*it2(end_column2) + 1; ...
                        it1(end_column2) + (n1+1)*it2(end_column2) + 2; ...
                        it1(end_column2) + (n1+1)*(it2(end_column2)+1) + 1; ...
                        ones(1,length(end_column2))];

    % Upper triangles |/ and \|
    t(1:4,(nt/2+1):2:nt)=[it1(1:2:(nt/2)) + (n1+1)*it2(1:2:(nt/2)) + 1; ... 
                          it1(1:2:(nt/2)) + (n1+1)*(it2(1:2:(nt/2))+1) + 2; ...
                          it1(1:2:(nt/2)) + (n1+1)*(it2(1:2:(nt/2))+1) + 1; ...
                          ones(1,length((nt/2+1):2:nt))];

    t(1:4,end_column1)=[it1(end_column2) + (n1+1)*it2(end_column2) + 2; ... 
                        it1(end_column2) + (n1+1)*(it2(end_column2)+1) + 2; ...
                        it1(end_column2) + (n1+1)*(it2(end_column2)+1) + 1; ...
                        ones(1,length(end_column2))];
end

% Shuffle t to get around qsort() bug on some platforms
r=rand(1,nt);[r,i]=sort(r);t=t(:,i);

