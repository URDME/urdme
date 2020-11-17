
%%
close all;

model = createpde(1);
% gd = [3 4 -1 1 1 -1 -1 -1 1 1]';
% sf = 'SQ1';
% ns = char(sf)';
% G = decsg(gd,sf,ns);
% [X,Y] = meshgrid(-1:2/120:1);
% X = X(:);
% Y = Y(:);
% nodes = [X';Y'];
% K = delaunayTriangulation(nodes');
% elements = K.ConnectivityList';
Nvoxels = 121;
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[L,M] = assema(P,T,1,1,0);
TR = triangulation(T(1:3,:)',P');
geometryFromMesh(model,TR.Points',TR.ConnectivityList');
specifyCoefficients(model,'m',0,'d',0,'c',1,'a',1,'f',0);
% figure;
% pdemPesh(model)
% axis equal
% state.time = 0;
FEM = assembleFEMatrices(model);
% applyBoundaryCondition(model,'edge',1:model.Geometry.NumEdges,'u',0);
% figure;
% spy(FEM.A)
% figure;
% spy(M);
thirdmatrix = round(M-FEM.A, 14);
figure;
title('Diff M')
spy(thirdmatrix);

% figure;
% spy(FEM.K)
% figure;
% spy(L);
thirdmatrix = round(L-FEM.K, 14);
figure;
title('Diff L')
spy(thirdmatrix);

% the (lumped) mass matrix gives the element volume
dM2 = full(sum(FEM.A,2));
ndofs = size(dM2,1);

% explicitly invert the lumped mass matrix and filter the diffusion matrix
[i,j,s] = find(FEM.K);
s = s./dM2(i);
%keep = find(s < 0); % (possibly removes negative off-diagonal elements)
keep = find(i ~= j); % (removes only the diagonal)
i = i(keep); j = j(keep); s = s(keep);

% rebuild L, ensuring that the diagonal equals minus the sum of the
% off-diagonal elements
L2 = sparse(i,j,s,ndofs,ndofs);
L2 = L2+sparse(1:ndofs,1:ndofs,-full(sum(L2,2)));

%% Find points that make up the boundary

[row,col] = find(ismember(T(1:3,:),idof1));

% indices to unique values in col
[~, ind] = unique(col, 'rows');
% duplicate indices
duplicate_ind = setdiff(1:size(col, 1), ind);
% % duplicate values
% duplicate_value = col(duplicate_ind);

cont_int_points = zeros(2,length(duplicate_ind));

for i = 1:length(duplicate_ind)
    cont_int_points(1,i) = T(row(duplicate_ind(i)-1),col(duplicate_ind(i)-1));
    cont_int_points(2,i) = T(row(duplicate_ind(i)),col(duplicate_ind(i)));
end
cont_int_points = sort(cont_int_points,1);
[~,inds] = unique(cont_int_points(1,:));
cont_int_points = cont_int_points(:,inds);

plot(P(1,cont_int_points),P(2,cont_int_points),'*')

%%

keep2 = find(ismember(Adof(jj_),idof1));
iii = reshape(ii(keep2),[],1); jjj_ = reshape(jj_(keep2),[],1);
idof1_moves = sort(max(Pr(sdof_m_(iii))-Pr(jjj_),0).*Drate_(2*VU(Adof(jjj_))+abs(U(Adof(jjj_)))+1))
m = mean(max(Pr(sdof_m_(ii))-Pr(jj_),0).*Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1))
maxmax = max(max(Pr(sdof_m_(ii))-Pr(jj_),0).*Drate_(2*VU(Adof(jj_))+abs(U(Adof(jj_)))+1))

%%
[P,E,T,gradquotient] = basic_mesh(1,121);
R = RobinMassMatrix2D(P,E);
figure; spy(R);

%%
for i = 1:size(Mgamma,1)
    nnz(Mgamma(i,:))
end

%%

