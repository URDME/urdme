close all;

% Problem parameters
Nvoxels = 41;
h = 2/(Nvoxels-1);
[x,y] = meshgrid(linspace(-1,1,Nvoxels));
x_T = x'; y_T = y';
f = @(x,y) (cos(2*pi*x).*sin(2*pi*y));
G_D = sparse(Nvoxels,Nvoxels);
rhs_pos_D = [round((Nvoxels-2)^2/2)];
rhs_pos = [round(Nvoxels^2/2)];
rhs_val = 1;
% figure(69); 
% surf(x,y,f(x,y))
% title('f = cos(2*pi*x).*sin(2*pi*y)');
subf = false;

% Make matrixes for Dirichlet
m_D = Nvoxels-2;
B_D = spdiags(bsxfun(@times,ones(m_D,1),[-1 2 -1]),[-1 0 1],m_D,m_D);
eyMat_D = speye(Nvoxels-2);
A_D = kron(B_D,eyMat_D) + kron(eyMat_D,B_D);
A_D = A_D/h^2;


% Make matrixes for Neumann
m_N = Nvoxels;
B_N = spdiags(bsxfun(@times,ones(m_N,1),[-1 2 -1]),[-1 0 1],m_N,m_N);
B_N(1,2) = -2; B_N(end,end-1) = -2;
eyMat_N = speye(Nvoxels);
A_N = kron(B_N,eyMat_N) + kron(eyMat_N,B_N);
A_N = A_N/h^2;

% Create system for Dirichlet
lhs_D = A_D;
LHS = lhs_D;
% rhs_D = reshape(f(x(2:end-1,2:end-1),y(2:end-1,2:end-1)),[],1);
rhs_D = 0*ones(size(A_D,1),1); rhs_D(rhs_pos_D) = rhs_D(rhs_pos_D) + rhs_val;
RHS = rhs_D;
solution_D = (LHS \ RHS);
G_D(2:end-1,2:end-1) = reshape(solution_D,Nvoxels-2,Nvoxels-2);
if subf
    figure(42);
    subplot(1,2,1);
    surf(x,y,G_D);
    cond_D = rcond(full(LHS));
    title(sprintf('Dirichlet\n rcond = %d', cond_D));
    figure(43);
    subplot(1,2,1);
    hold on;
    [px,py] = gradient(G_D);
    contour(x,y,G_D)
    quiver(x,y,px,py)
    hold off;
    title(sprintf('Dirichlet\n rcond = %d', cond_D));
else
    figure('Name',"Surf_Dirichlet");
    surf(x,y,G_D);
    cond_D = rcond(full(LHS));
    title(sprintf('Dirichlet\n rcond = %d', cond_D));
    figure('Name',"Gradient_Dirichlet");
    hold on;
    [px,py] = gradient(G_D);
    contour(x,y,G_D)
    quiver(x,y,px,py,'AutoScale','on')
    hold off;
    title(sprintf('Dirichlet\n rcond = %d', cond_D));
end

% Create system for Neumann
lhs_N = A_N;
LHS = lhs_N;
% rhs_N = reshape(f(x,y),[],1);
rhs_N = 0*ones(size(A_N,1),1); rhs_N(rhs_pos) = rhs_N(rhs_pos) + rhs_val;
% RHS = rhs_N;
RHS = rhs_N - mean(rhs_N);
solution_N = (LHS \ RHS);
c = mean(solution_N);
solution_N = solution_N - c;
G_N = reshape(solution_N,Nvoxels,Nvoxels);
if subf
    figure(42);
    subplot(1,2,2);
    surf(x,y,G_N);
    cond_N = rcond(full(LHS));
    title(sprintf('Neumann\n rcond = %d', cond_N));
    figure(43);
    subplot(1,2,2);
    hold on;
    [px,py] = gradient(G_N);
    contour(x,y,G_N)
    quiver(x,y,px,py)
    hold off;
    title(sprintf('Neumann\n rcond = %d', cond_N));
else
    figure('Name',"Surf_Neumann");
    surf(x,y,G_D);
    hold on;
    freezeColors;
    G_N(3:end-2,3:end-2) = NaN;
    surf(x,y,G_N + 0.00011);
    colormap('autumn');
    cond_N = rcond(full(LHS));
%     title(sprintf('Neumann\n rcond = %d', cond_N));
    axis tight, axis off;
    figure('Name',"Gradient_Neumann");
    hold on;
    [px,py] = gradient(G_N);
    contour(x,y,G_N)
    quiver(x,y,px,py,'AutoScale','on')
    hold off;
    title(sprintf('Neumann\n rcond = %d', cond_N));
end