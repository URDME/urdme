function ok = spin_ssa(ix)
%SPIN_SSA Tests for SSA solver.

% S. Engblom 2024-05-08 (data_time, ldata_time, gdata_time)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2018-02-15 (Minor revision, seeds and Nreplicas)
% S. Engblom 2017-02-22

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4};
stests = {'Spin_ssa #1' 'Spin_ssa #2' 'Spin_ssa #3' 'Spin_ssa #4'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SPIN_SSA (URDME/SSA)', ftests(ix), stests(ix));
unix('rm mexssa*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Four cells with one reaction.
%   X+Y --> Z.

ok = 1;

% Create a simple model.
umod.D = sparse(12,12);

% Initial vector u0.
u0 = [3 3 0; 3 3 0; 3 3 0; 2 1 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 1 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin1','seed',1234, ...
             'report',0,'solver','ssa');

xx = umod.U(:,length(umod.tspan));
ok = ok && sum([xx(1) xx(4) xx(7) xx(10)]) == 1;
ok = ok && sum([xx(2) xx(5) xx(8) xx(11)]) == 0;
ok = ok && sum([xx(3) xx(6) xx(9) xx(12)]) == 10;
ok = ok && all(umod.U(:) >=0 );

% Inline propensities.
K = [1 0 0]';
I = [1 2 1]';

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod = urdme(vmod,'seed',1234,'compile',0,'solver','ssa');
ok = ok && all(umod.U(:) == vmod.U(:));

vmod.u0 = cat(3,u0,u0,u0);
vmod = urdme(vmod,'seed',[1234 1234 1234],'compile',0,'solver','ssa');
ok = ok && norm(umod.U-vmod.U(:,:,1),'fro') == 0;
ok = ok && norm(umod.U-vmod.U(:,:,2),'fro') == 0;
ok = ok && norm(umod.U-vmod.U(:,:,3),'fro') == 0;

%-------------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Four cells with three reactions.
%   X+X --> Z, Y+Y --> Z and Z --> @.

ok = 1;
 
% Create a simple model. 
umod.D = sparse(12,12);

% Initial vector u0. 
u0 = [10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0 = u0;
 
% Stochiometric matrix N in sparse format.
umod.N = sparse([-2 -0 1; 0 -2 1; 0 0 -1]');

% Dependency graph G in sparse format.
umod.G = sparse([0 0 0 1 0 0; ...
                 0 0 0 0 1 0; ...
                 0 0 0 1 1 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin3','seed',29,'solver','ssa');

xx = umod.U(:,end);
ok = ok && sum(xx([1 4 7 10])) == 0;
ok = ok && sum(xx([2 5 8 11])) == 2;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [4 0 0; 4 0 0; 0 1 0]';
I = [1 1 1; 2 2 1; 1 1 3]';
S = sparse(4,3);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'seed',29,'compile',0,'solver','ssa');
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 Five cells with five reactions and four species.
%   D --> A, A+A --> B, B+B --> C, C+C --> A, A+B+C --> D.
%   Cubic reaction turned off for sd == 2.

ok = 1;

% Create the model.
umod.D = sparse(20,20);

% Initial vector u0.
u0 = [30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([0 0 0 0 1 0 0 0 1;
                 0 0 0 0 0 1 0 0 1;
                 0 0 0 0 0 0 1 0 1;
                 0 0 0 0 0 0 0 1 1;
                 0 0 0 0 0 1 1 1 1]);

% Output times tspan.
umod.tspan = [0:10:100];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(5) = 2;
umod.sd = sd;

umod = urdme(umod,'propensities','src/spin6','seed',20120303,'solver','ssa');

% Inline propensities.
K = [0 1 0;
     2 0 0;
     2 0 0;
     2 0 0]';

I = [1 1 4;
     1 1 1;
     2 2 1;
     3 3 1]';
% D -> A, A+A -> B, B+B -> C, C+C -> A, A+B+C -> D.
S = sparse(zeros(1,4));
S(1,1) = 2;

% using that 4 reactions are inline, the cubic one is first in src/spin6inline.c
vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'propensities','src/spin6inline', ...
             'seed',20120303,'solver','ssa');
ok = ok && all(umod.U(:) == vmod.U(:));

% check that non-empty D throws:
umod.D(1,1) = 1;
try
  umod = urdme(umod,'propensities','src/spin6','seed',20120303,'solver','ssa');
  ok = 0;
catch
  ok = ok;
end

%-------------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 Test of gdata_time, ldata_time.

ok = 1;

umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata_time' 'Gamma' 'gdata'},'src/spin4SIR.c');
umod.gdata = 0.2;
umod.data_time = [0 10 20 30 50];
umod.gdata_time = [1.5 3 0 1 5]*umod.gdata;

umod.tspan = 0:1:100;
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solver','ssa','seed',19443);

% switch to scalar ldata_time/ldata-fields instead
vmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata'},'src/spin4SIR.c');
vmod.ldata = 0.2;
vmod.data_time = [0 10 20 30 50];
vmod.ldata_time = permute([1.5 3 0 1 5]*umod.gdata,[1 3 2]);
vmod.tspan = 0:1:100;
vmod.D = sparse(3,3);
vmod.vol = 1e4;
vmod.sd = 1;
vmod.u0 = [umod.vol-100 100 0]';
vmod = urdme(vmod,'solver','ssa','seed',19443);

ok = ok && all(umod.U == vmod.U,'all');

% only ldata_time fields, but in 3 cells
wmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata_time'},'src/spin4SIR.c');
wmod.ldata = [0.2 0.2 0.2];
wmod.data_time = [0 10 20 30 50];
wmod.ldata_time = permute(repmat([[1.5 3   0   1   5]*umod.gdata; ...
                    0.2 0.2 0.2 0.2 0.2],[1 1 3]),[1 3 2]);
wmod.tspan = 0:1:100;
wmod.D = sparse(3*3,3*3);
wmod.vol = [1e4 1e4 1e4];
wmod.sd = [1 1 1];
wmod.u0 = repmat([wmod.vol(1)-100 100 0]',1,3);;

wmod = urdme(wmod,'solver','ssa','seed',19443);
ok = ok && all(wmod.U(1:3,:) == umod.U,'all');

unix('rm src/spin4SIR.c');

% $$$ figure(1), clf, plot(umod.tspan,umod.U);
% $$$ xline(umod.data_time);
% $$$ legend('S','I','R');

%-------------------------------------------------------------------------------
