function r = RobinLoadVector2D(P,T)
np = size(P,2);
nt = size(T,2);
r = zeros(np,1);

inds = [1,2;2,3;3,1];
len = [0.0167, 0.0167, 0.0167];
tmp = ones(1,3);
for K = 1:nt
    loc2glb = T(1:3,K); 
    x = P(1,loc2glb); % x-coordinates of triangle nodes
    y = P(2,loc2glb); % y-coordinates of triangle nodes
    len = [0.0167, 0.0167, 0.0167];
    tmp = ones(1,3);
    for ii = 1:length(loc2glb)
%        len(1,ii) = sqrt((x(1+mod(ii-1,length(loc2glb)))-x(1+mod(ii,length(loc2glb))))^2 ...
%            +(y(1+mod(ii-1,length(loc2glb)))-y(1+mod(ii,length(loc2glb))))^2);
       if x(inds(ii,1)) ~= x(inds(ii,2)) && y(inds(ii,1)) ~= y(inds(ii,2))
           len(1,ii) = 0;
           tmp(1,ii) = 0;
       end
    end
    rE = ([1 1 1].*(tmp+len)./(4*3))';
    r(loc2glb) = r(loc2glb) + rE;
end