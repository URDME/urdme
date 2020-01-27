function umod = schnakenberg(umod)
%Model file for the Schnakenberg reaction

% S. Engblom 2019-11-27 (Revision, augment model using rparse)
% S. Engblom 2017-02-20 (Revision, URDME 1.3, Comsol 5)
% Y. Saygun 2015-06-01 (Computational Stochastic Morphogenesis)

% transitions and rates
r1 = '@ > k1*vol > U';
r2 = 'U > k2*U > @';
r3 = '@ > k3*vol > V';
r4 = 'U+U+V > (k4/(vol*vol))*(U*(U-1)*V) > U+U+U';

umod = rparse(umod,{r1 r2 r3 r4},{'U' 'V'}, ...
              {'k1' 0.1 'k2' 1 'k3' 0.9 'k4' 1}, ...
              'schnakenberg.c');

% random initial number of species
rng(1702);
umod.u0 = [randi(500,1,numel(umod.vol)); ...
           randi(50,1,numel(umod.vol))];

% simulation time interval
if ~isfield(umod,'tspan')
  umod.tspan = 0:1:100;
end
