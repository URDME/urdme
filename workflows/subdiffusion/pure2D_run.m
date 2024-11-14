% Subdiffusion 2D example.
%   Single species, pure subdiffusion in 2D.

% S. Engblom 2017-04-10
% S. Engblom 2017-02-20 (Revision)
% S. Engblom 2014-05-09

if ~exist('report','var')
  report = 1;
end

clear umod;

% (1) create the geometry, here composed of a circle
C1 = [1 0 0 1]';
gd = [C1];
sf = 'C1';
ns = char('C1')';
G = decsg(gd,sf,ns);

% (2) create the mesh
[P,E,T] = initmesh(G,'hmax',0.04);
Nmesh = size(P,2);

% (3) assemble the diffusion part
umod = pde2urdme(P,T,{0.01});
D = umod.D;
umod.sd = round(umod.sd); % (unused)

% internal states approximation from pre-simulations
if ~exist('phi_','var')
  phi_ = 0.4; % default value
end
[gamma,A,p] = subdata(phi_);
if ~exist('FAC','var')
  FAC = 1; % default value
end
A = FAC*A;
Nstates = numel(gamma);

% initial- and equilibrium effective diffusion
[~,igammamax] = max(gamma);
gamma0 = gamma(igammamax);
gammabar = p'*gamma;

% given the internal state transfer matrix A, the diffusion
% coefficients gamma, the scalar diffusion operator D, compute the
% effective subdiffusion matrix
umod.D = subdmatrix(A,gamma,umod.D);

% initial vector u0: put particles in the middle, and in the fastest
% diffusing state
N0 = 1e2;
R = P(1,:).^2+P(2,:).^2;
[~,imid] = min(R);
umod.u0 = zeros(Nstates,Nmesh);
umod.u0(igammamax,imid) = N0;

% unused
umod.N = sparse(Nstates,0);
umod.G = sparse(0,Nstates);
umod.vol = ones(1,Nmesh);
umod.sd = ones(1,Nmesh);

% output times
Tend = 100;
t = linspace(0,Tend,101);
umod.tspan = t;

% solve
umod = urdme(umod,'report',report,'solver','nsm');
vmod = urdme(umod,'report',report,'solver','uds','solverargs',{{'odesolv',@ode23}});

% pure diffusion versions
wmod = vmod;
wmod.u0 = sum(wmod.u0,1);
wmod.N = sparse(1,0);
wmod.G = sparse(0,1);
wmod.D = gamma0*D;
wmod1 = urdme(wmod);
wmod.D = gammabar*D;
wmod2 = urdme(wmod);

% visualize
U = reshape(umod.U,Nstates,Nmesh,[]);
U = squeeze(sum(U,1));
V = reshape(vmod.U,Nstates,Nmesh,[]);
V = squeeze(sum(V,1));
W1 = wmod1.U;
W2 = wmod2.U;

if exist('plotting_off','var') && plotting_off, return; end
figure(1), clf
msdU = R*U./sum(U);
msdV = R*V./sum(V);
msdW1 = R*W1./sum(W1,1);
msdW2 = R*W2./sum(W2,1);

%loglog(t,msdU,'b'); hold on
loglog(t,msdV,'r','linewidth',2); hold on
loglog(t,msdW1,'k--');
loglog(t,msdW2,'k--');

exact = 2*0.01*t(end)*[gammabar gamma0];
axis([t(1) t(end) 0 exact(1)*2]);
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 400 150]);
set(gca,'xtick',[1e0 1e1],'ytick',[1e-2 1e-1],'yaxislocation','right');
h = xlabel('time');
set(h,'fontsize',8)
h = ylabel('MSD');
set(h,'fontsize',8)

% select suitable time window
select = find(msdV <= 0.99*msdW1 & ...
              msdV >= 1/0.99*msdW2 & ...
              10^0.5/FAC <= t & t <= 10^1.5/FAC);
% (in between conditions + avoid boundary effects)

% slope alpha
coeff = polyfit(log(t(select)),log(msdV(select)),1);
coeff = coeff(1)
scale = 0.2*msdV(select(1))/t(select(1))^coeff;
tt = t(select);
loglog(tt,scale*tt.^coeff,'b','linewidth',2);

return;

% cut-and-paste code

% produce table
phi = 0.05:0.05:0.4;
alpha = zeros(size(phi));
i = 0;
for phi_ = phi
  i = i+1;
  pure2D_run
  alpha(i) = coeff;
end
% $$$ save data/phi_vs_alpha phi alpha
load data/phi_vs_alpha
arr2latex([phi; alpha],{'.2f' '.2f'}', ...
          'hline','off','colspec','lrrrrrrrr', ...
          'rowlabel',{'$\phi$' '$\alpha$'},'label', ...
          'tab:phivsalpha', ...
          'caption',['Estimated values of subdiffusive constant ' ...
                    '$\alpha$ for different values of $\phi$. These ' ...
                    'values were obtained from least square fitting ' ...
                    'within a window of time and are therefore subject ' ...
                    'to some degree of uncertainty.'])

return;

% print to file
phi_ = 0.35;
FAC = 1;
pure2D_run
 
figure(1),
%print -depsc msd2d1.eps
% $$$ coeff =
% $$$    0.691112323449392

phi_ = 0.35;
FAC = 4;
pure2D_run

figure(1),
%print -depsc msd2d2.eps
% $$$ coeff =
% $$$    0.705386742011130
