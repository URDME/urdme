% Lotka-Volterra 3D example.
%   Species Y1 (prey) and Y2 (predator) compete in a periodic
%   space. This example shows how one can set up a simple 3D
%   simulation by 1D-methods.
%
%   This model is inspired from a similar one in [1].
%
%   Reference:
%     [1] S. S. Andrews and D. Bray: "Stochastic simulation of
%     chemical reactions with spatial resolution and single molecule
%     detail", Phys. Biol. 1(3):137--151 (2004).

% S. Engblom 2018-07-03

%% (1) geometry and diffusion operator

% diffusion rate
D_const = 2;

% build the 3D diffusion operator from 1D constructs
Ldim = [200 200 20]; % geometry dimensions
h = 4; % mesh discretization parameter

D = [];
mesh = cell(1,3);
for dim = 1:3
  % 1D diffusion operator
  mesh{dim} = h/2:h:Ldim(dim)-h/2; % voxel midpoints
  Ncells = numel(mesh{dim});
  e = ones(Ncells,1);
  D1 = spdiags([e -2*e e],-1:1,Ncells,Ncells);

  % periodic boundary
  D1(1,end) = e(1);
  D1(end,1) = e(end);
  
  % recursive construct
  if isempty(D)
    D = D1;
  else
    D = kron(D1,speye(size(D)))+kron(speye(size(D1)),D);
  end
end
D = 1/h^2*D; % final scaling

Ncells = size(D,1);
umod.vol = repmat(h^3,1,Ncells);
umod.sd = ones(1,Ncells);

% two species with same diffusion rate
umod.D = D_const*kron(D,speye(2));

%% (2) reactions
VOL = prod(Ldim);
umod = rparse_inline(umod, ...
                     {'Y1 > c1 > Y1+Y1', ...
                    'Y1+Y2 > c2 > Y2+Y2', ...
                    'Y2 > c3 > @'}, ...
                     {'Y1' 'Y2'},{'c1' 0.2 'c2' 100 'c3' 0.2});

%% (3) run the model

% initial conditions: 1000 Y1 and 1000 Y2 randomly distributed
Mspecies = 2;
umod.u0 = full([ ...
    sparse(1,randi(Ncells,1,1000),1,1,Ncells); ...
    sparse(1,randi(Ncells,1,1000),1,1,Ncells)]);

% simulate
if ~exist('tspan','var')
  tspan = 0:1:200;
end
umod = urdme(umod,'tspan',tspan,'seed',20180703,'report',0);

%% (4) visualize (can be turned off)
if ~exist('plotting_off','var') || ~plotting_off
  figure(1), clf,
  sumY1 = sum(umod.U(1:2:end,:),1);
  sumY2 = sum(umod.U(2:2:end,:),1);
  plot(tspan,sumY1,'r',tspan,sumY2,'g');
  legend('Y1','Y2');

  figure(2), clf
  plot(sumY1,sumY2,'k');
  xlabel('Y1'); ylabel('Y2');

  figure(3), clf,
  U = reshape(umod.U(:,end),2, ...
              numel(mesh{1}),numel(mesh{2}),numel(mesh{3}));
  U = sum(U,4); % sum along z-dimension
  Y1 = squeeze(U(1,:,:));
  Y2 = squeeze(U(2,:,:));
  [i,j,~] = find(Y1);
  % distribute individuals within voxel in a uniform random way:
  x = mesh{1}(i)+h*(rand(1,numel(i))-0.5);
  y = mesh{2}(j)+h*(rand(1,numel(j))-0.5);
  plot(x,y,'r.'); hold on,
  [i,j,~] = find(Y2);
  x = mesh{1}(i)+h*(rand(1,numel(i))-0.5);
  y = mesh{2}(j)+h*(rand(1,numel(j))-0.5);
  plot(x,y,'go');
  axis([0 Ldim(1) 0 Ldim(2)]);
  
  % animation
  figure(4), clf
  for k = 1:numel(umod.tspan)
    U = reshape(umod.U(:,k),2, ...
                numel(mesh{1}),numel(mesh{2}),numel(mesh{3}));
    U = sum(U,4); % sum along z-dimension
    Y1 = squeeze(U(1,:,:));
    Y2 = squeeze(U(2,:,:));
    [i,j,~] = find(Y1);
    % distribute individuals within voxel in a uniform random way:
    x = mesh{1}(i)+h*(rand(1,numel(i))-0.5);
    y = mesh{2}(j)+h*(rand(1,numel(j))-0.5);
    plot(x,y,'r.'); hold on,
    [i,j,~] = find(Y2);
    x = mesh{1}(i)+h*(rand(1,numel(i))-0.5);
    y = mesh{2}(j)+h*(rand(1,numel(j))-0.5);
    plot(x,y,'go'); hold off,
    axis([0 Ldim(1) 0 Ldim(2)]);
    drawnow;
  end
end
