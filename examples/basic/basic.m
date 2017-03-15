function umod = basic(umod)
%BASIC Model file for basic URDME example.

% S. Engblom 2013-03-27 (Minor revision)
% P. Bauer 2013-03-20

% sizes of problem
Mspecies = 3; % ordering is [X Y Z]
Ncells = numel(umod.vol); 

% stochiometric matrix N in sparse format
umod.N = sparse([-1  -1  1; ... % X+Y --> Z     
                  1   1 -1]');  % Z --> X+Y

% dependency graph G in sparse format
umod.G = sparse([1 1 0 1 1; ...
                 0 0 1 1 1]);
% 3 first columns correspond to diffusion events,
% 2 following columns to reactions...
% ...and the 2 rows are the two reactions

% initial number of species 
umod.u0 = zeros(Mspecies,Ncells); 

% distribute 50 X and 50 Y randomly
rng(123); % to get reproducible results
for i = 1:50
  cell = randi(Ncells,1,1);
  umod.u0(1,cell) = umod.u0(1,cell)+1;

  cell = randi(Ncells,1,1);
  umod.u0(2,cell) = umod.u0(2,cell)+1;
end

% output solution at times tspan
umod.tspan = 0:.1:100;
