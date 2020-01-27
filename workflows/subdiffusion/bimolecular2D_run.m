% Subdiffusion example with reactions.

% S. Engblom 2018-02-16 (minor revision)
% S. Engblom 2017-05-30 (added Collins-Kimball)
% S. Engblom 2017-04-10 (Major revision)
% S. Engblom 2017-02-20 (Revision)
% S. Engblom 2014-05-09

clear umod;

if ~exist('Kcase','var'), error('Define Kcase first.'); end

%% (1) create the geometry, here composed of a circle
C1 = [1 0 0 1]';
gd = [C1];
sf = 'C1';
ns = char('C1')';
G = decsg(gd,sf,ns);

%% (2) create the mesh
[P,E,T] = initmesh(G,'hmax',0.04);
Nmesh = size(P,2);

%% (3) assemble the diffusion part
umod = pde2urdme(P,T,{0.01 0.01 0.01});
umod.sd = round(umod.sd); % (unused)
vmod = umod;

% internal states approximation from pre-simulations
[gamma,A,p] = subdata(0.2);
Nstates = numel(gamma);

% initial- and equilibrium effective diffusion
[~,igammamax] = max(gamma);
gamma0 = gamma(igammamax);
gammabar = p'*gamma;

% given the internal state transfer matrix A, the diffusion
% coefficients gamma, the scalar diffusion operator D, compute the
% effective subdiffusion matrix
umod.D = subdmatrix(A,gamma/gammabar,umod.D,[1 1 0]);
% note: only the first two species (A and B) are expanded into
% internal states; also, the diffusion operators are scaled with
% gammabar to agree in steady state with the normal diffusion operator
% of vmod
ndofs = size(umod.D,1);

% initial vector u0: put particles uniformly, and in the fastest
% diffusing state
N0 = 1e5;
vmod.u0 = zeros(3,Nmesh);
vmod.u0(1,:) = full(sparse(1,randi(Nmesh,[1 N0]),1));
vmod.u0(2,:) = full(sparse(1,randi(Nmesh,[1 N0]),1));
% internal states version:
umod.u0 = zeros(Nstates+Nstates+1,Nmesh);
umod.u0([igammamax Nstates+igammamax],:) = vmod.u0(1:2,:);

%% (4) form the reaction model: a single bimolecular reaction
k0 = 1e-4;
[K,I,N,G] = rparse_inline({'A+B > k0 > C'},{'A' 'B' 'C'},{'k0' k0});
vmod.solverargs = {'K' K 'I' I};
vmod.N = N;
vmod.G = G;

% subdiffusive model for reactions Ai+Bj --> C with rate K0(i,j)
K0 = zeros(Nstates);   % K0(i,j) is the rate for Ai+Bj
switch Kcase
 case 1, K0(1,1) = 1;
 case 2, K0(Nstates,Nstates) = 1;
 case 3, K0 = (1:10)'*(1:10);
 case 4, K0 = (10:-1:1)'*(10:-1:1);
 case 5,
  % Collins-Kimball
  kb = k0;
  rhor = 0.2*k0;
  Dgamma = repmat(gamma/gammabar,1,Nstates);
  Dgamma = Dgamma+Dgamma';
  K0 = (kb*4*pi*Dgamma*rhor)./(kb+4*pi*Dgamma*rhor);
end
% normalization to be comparable to non-subdiffusive model
K0 = k0/(p'*K0*p)*K0;

[K,I,N] = rparse_inline({'A$i+B$j > K0_$i_$j > C'}, ...
                        {'A$i' 'B$j' 'C'}, ...
                        {'K0_$i_$j' K0}, ...
                        {'i' 1:Nstates 'j' 1:Nstates});

% remove all void reactions
ikeep = find(K(1,:));
K = K(:,ikeep);
I = I(:,ikeep);
N = N(:,ikeep);
umod.N = N;
umod.solverargs = {'K' K 'I' I};

% dependency graph G
umod.G = G_inline(K,I,umod.N);

% output times
Tend = 20;
t = linspace(0,Tend,101);
umod.tspan = t;
vmod.tspan = t;

% solve
umod = urdme(umod,'report',1);
if Kcase == 1
  vmod = urdme(vmod,'report',1);
end

% visualize
U = reshape(umod.U,[],Nmesh,numel(umod.tspan));
A = permute(sum(sum(U(1:Nstates,:,:),1),2),[3 1 2]);
B = permute(sum(sum(U(Nstates+1:2*Nstates,:,:),1),2),[3 1 2]);
C = permute(sum(U(end,:,:,:),2),[3 1 2]);
if Kcase == 1
  U0 = reshape(vmod.U,3,Nmesh,[]);
  C0 = permute(sum(U0(3,:,:),2),[2 3 1]);
end

return;

% cut-and-paste code:
Csave = cell(1,5);
for Kcase = 1:5
  bimolecular2D_run
  Csave{Kcase} = C;
end
save data/bimolecular2D t Csave C0

return;

% visualize the results using saved data:
load data/bimolecular2D

figure(1), clf,
loglog(t,Csave{1},'b.-',t,C0,'k--');
hold on;
for Kcase = 2:4
  switch Kcase
   case 2, loglog(t,Csave{Kcase},'b');
   case 3, loglog(t,Csave{Kcase},'r');
   case 4, loglog(t,Csave{Kcase},'r.-');
   case 5, loglog(t,Csave{Kcase},'g-.'); % very similar
  end
end

return;

% possibly print to file:
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 380 200]);
axis([t([1 end]) 1e0 1e6])
h = xlabel('time');
set(h,'fontsize',8)
h = ylabel('#Reactions');
set(h,'fontsize',8)

% print to file
figure(1),
print -depsc annihilation.eps
