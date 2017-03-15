function umod = brusselator(umod)
%Model file for the Brusselator reaction

% S. Engblom 2017-02-20 (Revision, URDME 1.3, Comsol 5)
% Y. Saygun - Stochastic Morphogenesis

% transitions and rates
r1 = '@ > k1*a*vol > U';
r2 = 'U > k2*b*U > V';
r3 = 'U+U+V > (k3/(vol*vol))*(U*(U-1)*V) > U+U+U';
r4 = 'U > k4*U > @';

[~,umod.N,umod.G] = rparse({r1 r2 r3 r4},{'U' 'V'}, ...
                           {'k1' 1 'k2' 1 'k3' 1 'k4' 1 ...
                    'a' 4.5 'b' 7.5}, ...
                           'brusselator.c');

% zero initial number of species
umod.u0 = zeros(size(umod.N,1),numel(umod.vol));

% simulation time interval
if ~isfield(umod,'tspan')
  umod.tspan = 0:10:100;
end
