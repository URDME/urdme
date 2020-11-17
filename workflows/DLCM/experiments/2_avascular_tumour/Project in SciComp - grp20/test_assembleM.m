[r,c] = find(L);

MM = sparse(size(L));
tic
for i = r'
    for j = c'
        x = P(1,[i j]); % x-coordinates of triangle nodes
        y = P(2,[i j]); % y-coordinates of triangle nodes
        len = sqrt((x(1)-x(2))^2 +(y(1)-y(2))^2);
        MM(i,j) = len*((i == j)*2 + 1)/6;
    end
end
toc