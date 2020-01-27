% Notch-delta-reporter model: RDME model.
%
%   The intent with this experiment is to match the 0-dimensional
%   model with an spatial stochastic model.

% S. Engblom 2017-09-12 (RDME)
% S. Engblom 2017-09-11 (Revision)
% S. Engblom 2016-12-05

% simulation interval
if ~exist('Tend','var')
  warning('Non-existing variable (Tend). Default used.');
  Tend = 40;
end

% solution recorded at this grid:
tspan = linspace(0,Tend,51);

Nvoxels = size(U,1);
if ~exist('sol_restart','var') || ~sol_restart
  % delta-notch-reporter DOFs
  Notch = spreplace(U,abs(1e-3*param.betaN+1e-4*param.betaN*randn(nnz(U),1)));
  Delta = spreplace(U,abs(1e-3*param.betaD+1e-4*param.betaD*randn(nnz(U),1)));
  Report = sparse(size(U,1),size(U,2));
else
  % alternative, start from previous:
  clear sol_restart % force explicit reuse
  Notch = Nsave{end};
  Delta = Dsave{end};
  Report = Rsave{end};
end

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Nsave = cell(1,numel(tspan));
Dsave = cell(1,numel(tspan));
Rsave = cell(1,numel(tspan));

% RDME-layer description
[P_,E_,T_] = RDMElayer;

% assemble the diffusion part

% setup URDME model in variables [N D R N_a N_b D_a D_b]
clear umod
Diff = 1;
umod = pde2urdme(P_,T_,{Diff Diff Diff 0 0 0 0});
% (note: variables [N_a N_b D_a D_b] are signal variables and are in a
% sense 0-dimensional)
% *** implementation possibility: use ldata[0,1,2,3] instead

% subdomains
umod.sd = floor(umod.sd);
Nvoxels_ = numel(umod.sd);

% scale so that each cell is unit total volume
umod.vol = umod.vol/sum(umod.vol);

% the scaling variables for (N,D) and R enables a comparision to the
% deterministic model
if ~exist('VolND','var') || ~exist('VolR','var')
  error('Define volume variables.');
end

% kinetic model
r = cell(1,9);
% (1) birth N, (2) degradation N (trans), (3) mutual degradation
% (N,D), (4) degradation N
r{1} = '@ > betaN*vol > N';
r{2} = 'N > N*(wa*D_a+wb*D_b)/kt > @';
r{3} = 'N+D > N*D/(kc*vol) > @';
r{4} = 'N > N > @';

% (5) birth D, downregulated by R, (6) degradation D (trans), (7)
% degradation D
r{5} = '@ > betaD*vol/(1.0+pow(R/(VolR*vol),2.0)) > D';
r{6} = 'D > D*(wa*N_a+wb*N_b)/kt > @';
r{7} = 'D > D > @';

% (8) birth R, upregulated by (trans-)D and N, (9) death R
r{8} = '@ > betaR*vol/(kRS/pow((qa*D_a+qb*D_b)/VolND*(N/(vol*VolND)),2.0)+1.0) > R';
r{9} = 'R > R > @';

spec = {'N' 'D' 'R' 'N_a' 'N_b' 'D_a' 'D_b'};
rate = {'VolND' VolND 'VolR' VolR ...
        'betaN' VolND*param.betaN 'betaD' VolND*param.betaD ...
        'betaR' VolR*param.betaR ...
        'kt' VolND*param.kt 'kc' VolND*param.kc 'kRS' param.kRS ...
        'wa' param.wa 'wb' param.wb 'qa' param.qa 'qb' param.qb};
[~,umod.N,umod.G] = rparse(r,spec,rate,'source/NDR_RDME.c');
umod.solver = 'nsm';
umod.propensities = 'source/NDR_RDME.c';
umod.seed = 1809;
umod.report = 0; % ***

% initial data
Usave(1) = {U};
Nsave(1) = {Notch};
Dsave(1) = {Delta};
Rsave(1) = {Report};

% distribute randomly inside each cell
usum = [round(VolND*full(Notch)'); ...
        round(VolND*full(Delta)'); ...
        round(VolR*full(Report)')];

umod.u0 = zeros(7,Nvoxels_,Nvoxels);
for i = 1:Nvoxels
  ind = floor(Nvoxels_*rand(1,usum(1,i)))+1;
  umod.u0(1,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
  ind = floor(Nvoxels_*rand(1,usum(2,i)))+1;
  umod.u0(2,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
  ind = floor(Nvoxels_*rand(1,usum(3,i)))+1;
  umod.u0(3,:,i) = full(sparse(1,ind,1,1,Nvoxels_));
end

for j = 2:numel(tspan)
  j

  % determine a suitable micro time-step
  dT = tspan(j)-tspan(j-1);
  X = full([Notch; Delta; Report]);
  dX = NDRrhs(0,X,Nj,Np,param);
  % allow at most 5% Euler step increase in a norm-wise sense:
  dtau = 0.05*norm(X)./norm(dX);
  % time-step is min of those, but using an integer number of steps:
  dt = min(dtau,dT);
  K = min(ceil(dT/dt),10); % *** to think of
  dt = dT/K;
  umod.tspan = [0 dt];

  % take K micro-steps
  for k = 1:K
    % update junctional and protrusional singalling contacts
    umod.u0(4,:,:) = reshape(repmat(usum(1,:)*Nj'-usum(1,:),Nvoxels_,1),1,Nvoxels_,[]);
    umod.u0(5,:,:) = reshape(repmat(usum(1,:)*Np'-usum(1,:),Nvoxels_,1),1,Nvoxels_,[]);
    umod.u0(6,:,:) = reshape(repmat(usum(2,:)*Nj'-usum(2,:),Nvoxels_,1),1,Nvoxels_,[]);
    umod.u0(7,:,:) = reshape(repmat(usum(2,:)*Np'-usum(2,:),Nvoxels_,1),1,Nvoxels_,[]);
    % (in reality this belongs to the DLCM-layer: somewhat an abuse of the
    % RDME-layer)

    % parse/compile/run NSM-solver:
    umod.seed = umod.seed+1;
    umod = urdme(umod);
    % will not parse/compile again:
    umod.parse = 0;
    umod.compile = 0;
 
    % prepare for next step
    umod.u0 = reshape(umod.U(:,end,:),[],Nvoxels_,Nvoxels);
    usum = permute(sum(umod.u0,2),[1 3 2]);
    % (note: umod.u0(4:7,:) is next overwritten)
  end

  % solution at grid point is now available
  Notch = sparse(usum(1,:)/VolND)';
  Delta = sparse(usum(2,:)/VolND)';
  Report = sparse(usum(3,:)/VolR)';
  Usave(j) = {U};
  Nsave(j) = {Notch};
  Dsave(j) = {Delta};
  Rsave(j) = {Report};
end

% display final frame
figure(3), clf,
i = numel(Usave);
NDRplot(Dsave{i},Usave{i},V,R);

% another attempt
mx = full(max(Dsave{end}(:)));
figure(8),
i = numel(Usave);
clf,
patch('Faces',R,'Vertices',V, ...
      'FaceColor',[0.9 0.9 0.9],'EdgeColor','none');
hold on,
patch('Faces',R,'Vertices',V,'EdgeColor',[0.7 0.7 0.7], ...
      'FaceVertexCData',log2(1/1024+full(Dsave{i}(:,1)/mx)),'FaceColor','flat');
colormap('bone')
colorbar;
axis equal, axis([-1 1 -1 1]);
set(gca,'xtick',[],'ytick',[]);
set(gca,'visible','off');
drawnow;

return;

% movie visualization
mx = full(max(Dsave{end}(:)));
figure(3),
M = struct('cdata',{},'colormap',{});
for i = 1:numel(Usave)
  clf,
  NDRplot(Dsave{i},Usave{i},V,R,mx);
  M(i) = getframe(gcf);
end

movie2gif(M,{M([1:2 end]).cdata}, ...
          'visual/NDR_RDME.gif', ...
          'delaytime',0.1,'loopcount',0);

return;

% look inside voxels
figure(3), h = gca;
while 1
  % fetch voxel
  axes(h);
  [x,y,button] = ginput(1);
  if button == 3, break; end
  [~,ii] = min((P(1,:)-x).^2+(P(2,:)-y).^2);
  ii

  % plot the RDME-layer content of that voxel
  u = reshape(umod.U(:,end),[],Nvoxels_,Nvoxels);
  u = u(:,:,ii);
  u(2,:) % Delta

  names = {'N' 'D' 'R'};
  for kk = 1:3
    figure(4+kk), clf,
    pdemesh(P_,E_,T_,u(kk,:)');
    title(names{kk});
  end
end

return;