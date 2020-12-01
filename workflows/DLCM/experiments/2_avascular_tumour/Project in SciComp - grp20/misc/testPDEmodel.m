
%%
close all;

model = createpde(1);
Nvoxels = 121;
% gd = [3 4 -1 1 1 -1 -1 -1 1 1]';
% sf = 'SQ1';
% ns = char(sf)';
% G = decsg(gd,sf,ns);
[P,E,T,grad] = flipped_mesh(Nvoxels);
pdemesh(P,E,T)
axis equal

% [P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
% [L,M] = assema(P,T,1,1,0);
% TR = triangulation(T(1:3,:)',P');
% geometryFromMesh(model,TR.Points',TR.ConnectivityList');
% specifyCoefficients(model,'m',0,'d',0,'c',1,'a',1,'f',0);
figure;
pdemesh(model)
axis equal
figure;
pdegplot(model.Geometry,'EdgeLabels','on')
axis equal
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
Nvoxels = 3;
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);
[L,dM,N,M] = dt_operators(P,T);
[L_orig,M_orig] = assema(P,T,1,1,0);
[V,R] = mesh2dual(P,E,T,'voronoi');
Mgamma_test = assemble_Mgamma(P,T);

alpha = [1e-3, 1e-1, 5e-1, 1e+0, 5e+0, 1e+1, 1e+3, 1e+4, 1e+5];
alpha_inv = 1./alpha;
adof = 1:floor(Nvoxels^2/2);
tdof = adof(end)+1:size(L,1);
Lai = fsparse(tdof,tdof,1,size(L));
Mgamma_test2_dM = Mgamma_test./dM;
diag_Mgamma = spdiags(Mgamma_test,0);
alpha_mat = fsparse(tdof,tdof,diag_Mgamma(tdof)',size(L));
i = 1;
for a_inv = alpha_inv

    % tic
    lhs_1 = (M \ (L_orig - Lai*L_orig + a_inv*Lai*Mgamma_test));
%     rhs_1 = full(fsparse([adof'; tdof'],1, ...
%         [ones(size(adof')) + 1./dM(adof); ones(size(tdof'))],[size(L,1) 1]));
    rhs_1 = ones(size(L,1),1) + full(fsparse(adof',1,1./dM(adof),[size(L,1) 1]));
    X_1 = lhs_1 \ rhs_1;
    % toc
    lhs_2 = L - Lai*L;
%     figure;spy(lhs);
    lhs_2 = lhs_2 + a_inv*Lai*Mgamma_test2_dM;
%     lhs = lhs + (Lai*Mgamma_test2 - alpha_mat + a_inv*alpha_mat);
%     figure;imagesc(lhs);
%     rhs_2 = full(fsparse([adof'; tdof'],1, ...
%         [ones(size(adof')) + 1./dM(adof); ones(size(tdof'))],[size(L,1) 1]));
    rhs_2 = ones(size(L,1),1) + full(fsparse(adof',1,1./dM(adof),[size(L,1) 1]));
%     figure;imagesc(rhs);
    X_2 = lhs_2 \ rhs_2;

%     figure;imagesc(lhs_1);
%     figure;imagesc(lhs_2);
%     figure;plot(X_1);
%     figure;plot(X_2);

    figure(1);
    subplot(3,3,i);   
    plot(adof, X_1(adof),'.-') %, 'DisplayName', sprintf('alpha = %d', 1/a_inv));
    hold on;
    plot(tdof, X_1(tdof),'.-')
    title(sprintf('alpha = %d', 1/a_inv));
    grid on;
    
    figure(2);
    subplot(3,3,i);   
    plot(adof, X_2(adof),'.-') %, 'DisplayName', sprintf('alpha = %d', 1/a_inv));
    hold on;
    plot(tdof, X_2(tdof),'.-')
    title(sprintf('alpha = %d', 1/a_inv));
    grid on;
    i = i + 1;
end
% legend;
% figure; spy(Mgamma_test1);
% figure; spy(Mgamma_test2);
% figure; spy(L);
% figure; spy(M);
% figure; pdemesh(P,E,T)
% 
% Mgamma_test1(adof(end),:)
% Mgamma_test2(adof(end),:)
% M(adof(end),:)
