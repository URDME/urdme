function ok = spin(ix)
%SPIN Tests for URDME.

% S. Engblom 2024-05-08 (data_time, ldata_time, gdata_time)
% S. Engblom 2022-04-13 (Revision, urdme_validate_model)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2019-11-06 (Revision, URDMEstate_t)
% S. Engblom 2018-02-10 (Revision, Nreplicas)
% S. Engblom 2017-03-03 (Revision, rparse_inline)
% S. Engblom 2017-03-01 (Revision, rparse)
% S. Engblom 2017-02-18 (Revision, inline propensities)
% S. Engblom 2017-02-15 (Minor revision)
% P. Bauer 2012-08-31
% J. Cullhed 2008-06-16

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_spin0 ...
          @l_spin1 @l_spin2 @l_spin3 @l_spin4 @l_spin5 @l_spin6 @l_spin7 ...
          @l_spin8 @l_spin9 @l_spin10 @l_spin11 @l_spin12 @l_spin13 ...
          @l_spin14};
stests = {'Spin #0' ...
          'Spin #1' 'Spin #2' 'Spin #3' 'Spin #4' 'Spin #5' 'Spin #6' ...
          'Spin #7' 'Spin #8' 'Spin #9' 'Spin #10' 'Spin #11' ...
          'Spin #12' 'Spin #13' 'Spin #14'};

if nargin == 0
  ix = 1:size(ftests,2);
else
  ix = ix+1; % zero-offset, this spin only
end
ok = runtest('SPIN (URDME)',ftests(ix),stests(ix));

%-------------------------------------------------------------------------------
function ok = l_spin0
%L_SPIN0 Test of the intended high-level use of URDME.

ok = 1;

% (1) build mesh and diffusion operator
h = 1;
mesh = h/2:h:10-h/2; % voxel midpoints

% 1D-diffusion with periodic boundaries
Ncells = numel(mesh);
e = ones(Ncells,1);
D = spdiags([e -2*e e],-1:1,Ncells,Ncells);
D(1,end) = e(1);
D(end,1) = e(end);
D_const = 0.1;
D = D_const/h^2*D;

umod.vol = repmat(h^3,1,Ncells);
umod.sd = ones(1,Ncells);

% three species with the same diffusion rate, 7 without diffusion
umod.D = kron(D,spdiags([1 1 1 0 0 0 0 0 0 0]',0,10,10));

% (2) build reactions (including reaction counters R_r)
% check that the ghost rate is noticed:
warning('off','rparse:ghost_rates');
umod = rparse(umod, ...
              {'@ > k1*vol > X + R1' ...
               '@ > k1*vol > Y + R2' ...
               'X > mu1*X > R3' ...
               'Y > mu2*Y > R4' ...
               'X+Y > kk*X*Y/vol > Z + R5' ...
               'Z > kki*Z > X+Y + R6' ...
               'Z > mu3*Z > R7'}, ...
              [{'X' 'Y' 'Z'} seqexpand('R$r',{'r' 1:7})], ...
              {'k1' 0.1 'k2' 0.2 'mu1' 'gdata' 'mu2' 'gdata' ...
               'kk' 0.001 'kki' 'ldata' 'mu3' 0.01}, ...
              'src/spin0.c');
[~,wid] = lastwarn;
warning('on','rparse:ghost_rates');
ok = ok && strcmpi(wid,'rparse:ghost_rates');

% insert reaction rates into gdata/ldata fields
ixmu1 = find(strcmp(umod.gdata,'mu1'));
ixmu2 = find(strcmp(umod.gdata,'mu2'));
mu = zeros(1,2);
mu([ixmu1 ixmu2]) = [0.01 0.02];
umod.gdata = mu;

ixkki = find(strcmp(umod.ldata,'kki'));
kki = zeros(1,1);
kki(ixkki) = [0.002];
umod.ldata = repmat(kki,1,Ncells);

% (3) initial conditions: 100 X and 100 Y randomly distributed
Mspecies = 10;
umod.u0 = full([ ...
    sparse(1,randi(Ncells,1,100),1,1,Ncells); ...
    sparse(1,randi(Ncells,1,100),1,1,Ncells); ...
    sparse(8,Ncells)]);

% (4) simulate
umod = urdme(umod,'seed',201911,'tspan', 0:10:100,'report',0);

% manually figure out what happened by using the counters R_r
sumU = squeeze(sum(reshape(umod.U,10,[],11),2));
manual = [sumU(:,1)+umod.N*sumU(4:end,:)];
ok = ok && norm(sumU-manual,'fro') == 0;

% cleanup
unix('rm src/spin0.c');
unix('rm mexnsm*.mex*');

% (5) check warning messages
try
  warn = 'URDME:noMake';
  warning('off',warn);
  umod = urdme(umod,'solver','nms');
  ok = 0;
  warning('on',warn);
catch
  warning('on',warn);
  [~,lastid] = lastwarn;
  ok = ok && strcmp(lastid,warn);
end
try
  warn = 'URDME:noMex';
  warning('off',warn);
  umod = urdme(umod,'solver','nms','compile',false);
  ok = 0;
  warning('on',warn);
catch
  warning('on',warn);
  [~,lastid] = lastwarn;
  ok = ok && strcmp(lastid,warn);
end

%-------------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Test of reaction-diffusion in four cells with one reaction.
%   X+Y --> Z

ok = 1;

% Create a simple model.
B = [ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2) = 0;
B(7:9,4) = 0;

D = spdiags(B,[-6 -3 0 3 6],12,12);
umod.D = D';

% Initial vector u0.
u0 = [3 3 0; 3 3 0; 3 3 0; 2 1 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 1 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other inputs
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin1','seed',1234);

xx = umod.U(:,length(umod.tspan));
ok = ok && sum([xx(1) xx(4) xx(7) xx(10)]) == 1;
ok = ok && sum([xx(2) xx(5) xx(8) xx(11)]) == 0;
ok = ok && sum([xx(3) xx(6) xx(9) xx(12)]) == 10;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [1 0 0]';
I = [1 2 1]';

vmod = umod;
vmod.inline_propensities.K = K;
vmod.inline_propensities.I = I;
vmod.compile = 0; % don't re-compile
vmod = urdme(vmod);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Test of reaction-diffusion in four cells with two reactions.
%   X+X --> Z and Y+Y --> Z

ok = 1;

% create model
B = [ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2) = 0;
B(7:9,4) = 0;
D = spdiags(B,[-6 -3 0 3 6],12,12);
umod.D = D';

% initial vector u0
u0 = [10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0 = u0;

% stochiometric matrix N
umod.N = sparse([-2 -0 1; -0 -2 1]');

% dependency graph G
umod.G = sparse([1 0 0 1 0; 0 1 0 0 1]);

% output times
umod.tspan = [0:10 70:10:200];

% inputs cell, vol, pos and sd
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin2','seed',123);

xx = umod.U(:,length(umod.tspan));
ok = ok && sum([xx(1) xx(4) xx(7) xx(10)]) == 0;
ok = ok && sum([xx(2) xx(5) xx(8) xx(11)]) == 0;
ok = ok && sum([xx(3) xx(6) xx(9) xx(12)]) == 50;
ok = ok && all(umod.U(:) >= 0);

% inline propensities
K = [2 0 0; 2 0 0]'*2;
I = [1 1 1; 2 2 1]';
S = sparse(4,2);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',[]);
vmod = urdme(vmod,'seed',123,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));

% repeat with multiple replicas
vmod.u0 = cat(3,u0,u0,u0);
vmod = urdme(vmod,'seed',123,'compile',0);
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,1),[],1));
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,2),[],1));
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,3),[],1));

% multiple seeds
vmod = urdme(vmod,'seed',[123 123 123],'compile',0);
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,1),[],1));
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,2),[],1));
ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,3),[],1));

% somewhat more detail
vmod.u0 = cat(3,u0+1,u0+2,u0+3);
vmod = urdme(vmod,'seed',[123 123 123]+(1:3),'compile',0);
for i = 1:3
  umod.u0 = u0+i;
  umod = urdme(umod,'seed',123+i,'compile',0);
  ok = ok && all(umod.U(:) == reshape(vmod.U(:,:,i),[],1));
end
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 Test of reaction-diffusion in four cells with three reactions.
%   X+X --> Z, Y+Y --> Z and Y --> Z.

ok = 1;
 
% Create a simple model. 
B = [ones(12,2) -2*ones(12,1) ones(12,2)]; 
B(4:6,2) = 0; 
B(7:9,4) = 0;
D = spdiags(B,[-6 -3 0 3 6],12,12); 
umod.D = D';

% Initial vector u0. 
u0 = [10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0 = u0;
 
% Stochiometric matrix N in sparse format.
umod.N = sparse([-2 -0 1; 0 -2 1; 0 0 -1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 1 1 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin3','seed',29);

xx = umod.U(:,length(umod.tspan));
ok = ok && sum([xx(1) xx(4) xx(7) xx(10)]) == 0;
ok = ok && sum([xx(2) xx(5) xx(8) xx(11)]) == 0;
ok = ok && sum([xx(3) xx(6) xx(9) xx(12)]) == 0;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [4 0 0; 4 0 0; 0 1 0]';
I = [1 1 1; 2 2 1; 1 1 3]';
S = sparse(4,3);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'seed',29,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 Test of reaction-diffusion in four cells with three reactions.
%   X+X --> Z and Y+Y --> Z.

ok = 1;

% Create a simple model. 
B = [ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2) = 0;
B(7:9,4) = 0;
D = spdiags(B,[-6 -3 0 3 6],12,12);
D(3:3:end,:) = 0;
D(:,3:3:end) = 0;
umod.D = D';

% Initial vector u0. 
u0 = [10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-2 0 1; 0 -2 1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 0 0 1 0; 0 1 0 0 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(1) = 2;
umod.sd = sd;

umod = urdme(umod,'propensities','src/spin4','report',0,'seed',178);

xx = umod.U(:,length(umod.tspan));
ok = ok && xx(3) == 50;
ok = ok && all(umod.U(:) >= 0);

K = [2 0 0; 2 0 0]'*2;
I = [1 1 1; 2 2 1]';
S = sparse([0 0; 1 1]);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'seed',178,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 Four cells with one reaction and diffusion in all cells except one.
%   X+Y --> Z.

ok = 1;

% Create a simple model.
B = [ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2) = 0;
B(1:3,3) = 0;
B(4:9,4) = 0;
B(7:9,5) = 0;
D = spdiags(B,[-6 -3 0 3 6],12,12);
umod.D = D';

% Initial vector u0.
u0 = [3 3 0; 3 3 0; 3 3 0; 2 1 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 0 1]);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin5','report',0,'seed',175);

xx = umod.U(:,length(umod.tspan));
ok = ok && xx(1) == 1 && xx(3) == 10;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [1 0 0]';
I = [1 2 1]';

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod = urdme(vmod,'seed',175);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin6
%L_SPIN6 Five cells with five reactions and four species.
%   D --> A, A+A --> B, B+B --> C, C+C --> A, A+B+C --> D.
%   Cubic reaction turned off for sd == 2.

ok = 1;

% Create the model.
B = [ones(20,2) -3*ones(20,1) ones(20,4)];
B(9:12,1) = 0;
B(5:8,2) = 0;
B(13:16,2) = 0;
B(17:20,3) = 0;
B(9:12,4) = 0;
B(13:16,6) = 0;
D = spdiags(B,[-8 -4 0 4 8 12 16],20,20);
umod.D = D';

% Initial vector u0.
u0 = [30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([0 0 0 1 1 0 0 0 1;
                 1 0 0 0 1 1 0 1 1;
                 0 1 0 0 0 1 1 0 1;
                 0 0 1 0 0 0 1 1 1;
                 1 1 1 0 1 1 1 1 1]);

% Output times tspan.
umod.tspan = [0:0.1:100];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(5) = 2;
umod.sd = sd;

umod = urdme(umod,'propensities','src/spin6','seed',20120303);

ok = ok && sum(umod.U(1:16,end)) == 0;
ok = ok && sum(umod.U(17:20,end)) == 79;
ok = ok && all(umod.U(:) >= 0);

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

% using that 4 reactions are inline, the cubic one is first in
% src/spin6inline.c
vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'propensities','src/spin6inline', ...
             'seed',20120303);
ok = ok && all(umod.U(:) == vmod.U(:));

% build the same model from scratch using rparse/rparse_inline
clear umod;
umod = rparse([],{'A+B+C > A*B*C/vol/vol > D'}, ...
              {'A' 'B' 'C' 'D'},{},'src/spin6inline_rparse');
umod = rparse_inline(umod,{'D > one > A' 'A+A > two > B' ...
                    'B+B > two > C' 'C+C > two > A'}, ...
                     {'A' 'B' 'C' 'D'},{'one' 1.0 'two' 2.0});
umod.inline_propensities.S(1,1) = 2;

umod.D = D';
umod.u0 = u0;
umod.tspan = [0:0.1:100];
umod.vol = ones(1,size(u0,2));
umod.sd = sd;
umod = urdme(umod,'seed',20120303);
ok = ok && all(umod.U(:) == vmod.U(:));

% in the other order:
clear umod;
umod = rparse_inline([],{'D > one > A' 'A+A > two > B' ...
                    'B+B > two > C' 'C+C > two > A'}, ...
                     {'A' 'B' 'C' 'D'},{'one' 1.0 'two' 2.0});
umod.inline_propensities.S(1,1) = 2;
umod = rparse(umod,{'A+B+C > A*B*C/vol/vol > D'}, ...
              {'A' 'B' 'C' 'D'},{},'src/spin6inline_rparse');

umod.D = D';
umod.u0 = u0;
umod.tspan = [0:0.1:100];
umod.vol = ones(1,size(u0,2));
umod.sd = sd;
umod = urdme(umod,'seed',20120303);
ok = ok && all(umod.U(:) == vmod.U(:));

unix('rm mexnsm*.mex*'); % cleanup
unix('rm src/spin6inline_rparse.c'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin7
%L_SPIN7 Six cells with three reactions and three species.
%   X+Y --> Z, X --> Y, Y --> X.

ok = 1;

% Create the model.
B = [-1*ones(18,1) ones(18,1)];
B(16:18,1) = 0;
D = spdiags(B,[0 3],18,18);
umod.D = D';

% Initial vector u0.
u0 = [1000 1000 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 1; -1 1 0; 1 -1 0]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 0 1 1 1;
                 1 0 0 1 1 1;
                 0 1 0 1 1 1]);

% Output times tspan.
umod.tspan = [0:0.3:101];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin7','report',0,'seed',99);

xx = umod.U(:,length(umod.tspan));
ok = ok && xx(18) == 1000;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [1 0 0;
     0 1 0;
     0 1 0]';
I = [1 2 1;
     1 1 1;
     1 1 2]';

vmod = umod;
vmod.inline_propensities = struct('I',I,'K',K);
vmod = urdme(vmod,'seed',99,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin8
%L_SPIN8 Four cells with three reactions and three species.
%   Diffusion is in opposite directions for x,w and y,z.
%   W+X --> Z, 100*(W+Y) --> Z, X+Y --> Z.

ok = 1;

% Create the model.
B1 = [0 1 1 0 0 1 1 0 0 1 1 0 0 0 0 0]';
B2 = [-1 0 0 -1*ones(1,9) 0 -1 -1 0]';
B3 = [0 0 0 0 1 0 0 1 1 0 0 1 1 0 0 1]';
B = [B1 B2 B3];
D = spdiags(B,[-4 0 4],16,16);
umod.D = D';

% Initial vector u0.
u0 = [100 0 0 100; 0 0 0 0; 0 0 0 0; 0 1000 100 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 0 1; 1 -1 -1 0; -1 1 -1 0]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 0 0 1 1 1;
                 0 1 1 0 1 1 1;
                 1 0 1 0 1 1 1]);

% Output times tspan.
umod.tspan = [0:0.321:101];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin8','report',0,'seed',899);

xx = umod.U(:,length(umod.tspan));
ok = ok && xx(2) == 700 && xx(end) == 300;
ok = ok && all(umod.U(:) >= 0);

K = [1 0 0;
     100 0 0;
     1 0 0]';
I = [1 2 1;
     2 3 1;
     1 3 1]';

vmod = umod;
vmod.inline_propensities = struct('I',I,'K',K);
vmod = urdme(vmod,'seed',899,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin9
%L_SPIN9 Four cells with one reaction and one species.
%   Y -> Z.

ok = 1;

% Create the model.
D = sparse([-2 1 1 0; 1 -2 0 1; 1 0 -2 1; 0 1 1 -2]);
umod.D = D';

% Initial vector u0.
u0 = [1; 100; 3; 5]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1]);

% Output times tspan.
umod.tspan = [0:0.321:101];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin9','report',0,'seed',1899);

xx = umod.U(:,length(umod.tspan));
ok = ok && any(xx == 0);
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [0 1 0]';
I = [1 1 1]';

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod = urdme(vmod,'seed',1899,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin10
%L_SPIN10 Test of random seed.
%   Compare two models with the same random seed and see if the result
%   is the same.
%   X+Y --> Z, X --> Y and Y --> Z.

ok = 1;

% Create the model.
B = [-1*ones(18,1) ones(18,1)];
B(16:18,1) = 0;
D = spdiags(B,[0 3],18,18);
umod.D = D';

% Initial vector u0.
u0 = [1000 1000 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([-1 -1 1; -1 1 0; 1 -1 0]');

% Dependency graph G in sparse format.
umod.G = sparse([1 1 0 1 1 1;
                 1 0 0 1 1 1;
                 0 1 0 1 1 1]);

% Output times tspan.
umod.tspan = [0:0.3:101];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));
umod.ldata = ones(0,size(u0,2));
umod.gdata = [];

umod1 = urdme(umod,'propensities','src/spin10','seed',20120329,'report',0);
umod2 = urdme(umod,'propensities','src/spin10','seed',20120329,'report',0);

ok = ok && all(all(umod1.U == umod2.U));
ok = ok && all(umod1.U(:) >= 0);
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin11
%L_SPIN11 Test of empty inline propensities.

% model from L_SPIN6
% D -> A, A+A -> B, B+B -> C, C+C -> A, A+B+C -> D.
ok = 1;

% create the model
B = [ones(20,2) -3*ones(20,1) ones(20,4)];
B(9:12,1) = 0;
B(5:8,2) = 0;
B(13:16,2) = 0;
B(17:20,3) = 0;
B(9:12,4) = 0;
B(13:16,6) = 0;
D = spdiags(B,[-8 -4 0 4 8 12 16],20,20);
umod.D = D';

u0 = [30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0 = u0;

umod.N = sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');
umod.G = sparse([0 0 0 1 1 0 0 0 1;
                 1 0 0 0 1 1 0 1 1;
                 0 1 0 0 0 1 1 0 1;
                 0 0 1 0 0 0 1 1 1;
                 1 1 1 0 1 1 1 1 1]);

umod.tspan = [0:0.1:100];

umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(5) = 2;
umod.sd = sd;
umod.ldata = ones(0,size(u0,2));
umod.gdata = [];

% various errors
try
  umod.inline_propensities = struct('K','I','I',zeros(3,0), ...
                                    'S',sparse(5,0));
  umod = urdme(umod,'propensities','src/spin6', ...
               'seed',20120303,'report',0);
  ok = 0;
catch
  ok = ok;
end
try
  umod.inline_propensities = struct('K',zeros(3,0),'I',zeros(3,0), ...
                                    'W',sparse(5,0));
  umod = urdme(umod,'propensities','src/spin6','compile',0, ...
               'seed',20120303,'report',0);
  ok = 0;
catch
  ok = ok;
end
try
  umod.inline_propensities = struct('K',zeros(3,0),'I',zeros(3,0), ...
                                    'S',sparse(5,1));
  umod = urdme(umod,'propensities','src/spin6','compile',0, ...
               'seed',20120303,'report',0);
  ok = 0;
catch
  ok = ok;
end

% final call, just empty inline propensities
umod.inline_propensities = ...
    struct('S',sparse(5,0),'K',zeros(3,0),'I',zeros(3,0));
umod = urdme(umod,'propensities','src/spin6','compile',1, ...
             'seed',20120303,'report',0);

ok = ok && sum(umod.U(:,end)) == 79;
ok = ok && all(umod.U(:) >= 0);
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin12
%L_SPIN12 Test of RPARSE, RPARSE_INLINE.

% model from L_SPIN6
% D -> A, A+A -> B, B+B -> C, C+C -> A, A+B+C -> D.
ok = 1;

umod = rparse([],{'D > sd != 2 ? D : 0.0> A','A+A > 1.0*A*(A-1)/vol > B', ...
                  'B+B>k*B*(B-1)/vol > C','C+C > 1*C*(C-1)/vol > A', ...
                  'A+B+C > k*A*B*C/(vol*vol)>D'},{'A' 'B' 'C' 'D'}, ...
                 {'k' 1}, ...
                 'src/spin6rparse.c');

% create the rest of the model
B = [ones(20,2) -3*ones(20,1) ones(20,4)];
B(9:12,1) = 0;
B(5:8,2) = 0;
B(13:16,2) = 0;
B(17:20,3) = 0;
B(9:12,4) = 0;
B(13:16,6) = 0;
D = spdiags(B,[-8 -4 0 4 8 12 16],20,20);
umod.D = D';
u0 = [30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0 = u0;

N0 = sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');
G0 = sparse([0 0 0 1 1 0 0 0 1;
             1 0 0 0 1 1 0 1 1;
             0 1 0 0 0 1 1 0 1;
             0 0 1 0 0 0 1 1 1;
             1 1 1 0 1 1 1 1 1]);

ok = ok && all(umod.N(:) == N0(:)) && all(umod.G(:) == G0(:));

umod.tspan = [0:0.1:100];

umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(5) = 2;
umod.sd = sd;
umod.ldata = ones(0,size(u0,2));
umod.gdata = [];

umod = urdme(umod,'propensities','src/spin6','seed',303);
vmod = urdme(umod,'propensities','src/spin6rparse');

ok = ok && sum(umod.U(:) == vmod.U(:));

% Inline propensities.
K0 = [0 1 0;
      2 0 0;
      2 0 0;
      2 0 0]';
I0 = [1 1 4;
      1 1 1;
      2 2 1;
      3 3 1]';

vmod = rparse_inline([],{ ...
    'D > k1 > A', ...
    'A+A > k2 > B', ...
    'B+B > k3 > C', ...
    'C+C > k4 > A'}, ...
                          {'A' 'B' 'C' 'D'}, ...
                          {'k1' 1 'k2' 2 'k3' 2 'k4' 2});

ok = ok && all(vmod.inline_propensities.K(:) == K0(:));
ok = ok && all(vmod.inline_propensities.I(:) == I0(:));
ok = ok && all(vmod.N(:) == reshape(umod.N(:,1:4),[],1));
ok = ok && all(vmod.G(:) == reshape(umod.G(1:4,1:end-1),[],1));

% cleanup
unix('rm src/spin6rparse.c');

% this was a bug previously:

% a species is a subset of another species

% check that the ghost species is noticed:
warning('off','rparse:ghost_species');
umod = rparse([],{'ABCD+C > ABCD*A*AB*C*AB*ABC*ABCD*AB > AB+ABC'}, ...
              {'C' 'ABC' 'ABCD' 'AB' 'A' 'B'},{},'~');
[~,wid] = lastwarn;
warning('on','rparse:ghost_species');
ok = ok && strcmpi(wid,'rparse:ghost_species');

F1 = umod.private.rp.code;
i = strfind(F1,'return'); i = i(1);
ok = ok && strcmp(F1(i:(i+97)), ...
                  'return xstate[ABCD]*xstate[A]*xstate[AB]*xstate[C]*xstate[AB]*xstate[ABC]*xstate[ABCD]*xstate[AB];');

% a species is a subset of a rate
% check that the ghost rate is noticed:
warning('off','rparse:ghost_rates');
umod = rparse([],{'A+BA+ABA > BB*A*BA*AA*ABA > B+BA'}, ...
              {'B' 'A' 'BA' 'ABA'},{'BB' 1 'AA' 2 'AB' 3},'~');
[~,wid] = lastwarn;
warning('on','rparse:ghost_rates');
ok = ok && strcmpi(wid,'rparse:ghost_rates');
F2 = umod.private.rp.code;
i = strfind(F2,'return'); i = i(1);
ok = ok && strcmp(F2(i:(i+45)),'return BB*xstate[A]*xstate[BA]*AA*xstate[ABA];');

% double use of '>'
umod = rparse([],{'X+Y > X > Y ? X*Y : 0.0 > XY'}, {'X' 'Y' 'XY'},{}, ...
              '~');
F3 = umod.private.rp.code;
i = strfind(F3,'return'); i = i(1);
ok = ok && strcmp(F3(i:(i+50)),'return xstate[X]>xstate[Y]?xstate[X]*xstate[Y]:0.0;');
unix('rm mexnsm*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok = l_spin13
%L_SPIN13 Test of URDME_VALIDATE_MODEL.

ok = 1;

umod = rparse([],{'X > k1 > @' ...
                  'X+Y > k2*X > @' ...
                  'X+X > k3*X > @' ...
                  'X+X > k4*(X-1) > @' ...
                  '@ > kk > Y'},{'X' 'Y'}, ...
              {'k1' 'ldata' ...
               'k2' 'ldata' ...
               'k3' 'ldata' ...
               'k4' 'ldata' ...
               'kk' 1},'src/spin13validate.c');
Ncells = 1;
umod.vol = ones(1,Ncells);
umod.tspan = [0 100];
umod.u0 = 10*ones(2,Ncells);
umod.D = sparse(2*Ncells,2*Ncells);
umod.sd = ones(1,Ncells);

I = eye(4);
for i = 1:4
  try
    umod.ldata = I(:,i);
    urdme_validate_model(umod);
    ok = 0;
  catch
    % check that it found rate/reaction #i:
    ok = ok && strfind(lasterr,sprintf('#%d',i));
  end
end
unix('rm src/spin13validate.c'); % cleanup
unix('rm *urdme_validate_model*mexrhs.mex*');
unix('rm *urdme_validate_model*mexjac.mex*');

%-------------------------------------------------------------------------------
function ok = l_spin14
%L_SPIN14 Test of RPARSE: gdata_time, ldata_time.

ok = 1;

% ldata
umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata' 'Gamma' 'ldata'},'src/spin14SIR.c');
umod.ldata = [1.5*0.2 0.2]';
umod.tspan = [0:1:100 1000];
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solve',1,'seed',1493);
U1 = umod.U;

% gdata
umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata' 'Gamma' 'gdata'},'src/spin14SIR.c');
umod.gdata = [1.5*0.2 0.2];
umod.tspan = [0:1:100 1000];
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solve',1,'seed',1493);
U2 = umod.U;
ok = ok && all(U1 == U2,'all');

% gdata_time
umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata_time' 'Gamma' 'gdata'},'src/spin14SIR.c');
umod.gdata = 0.2;
umod.data_time = [0 10 20 30 50];
umod.gdata_time = [1.5 1.5 1.5 1.5 1.5]*umod.gdata;
umod.tspan = [0:1:100 1000];
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solve',1,'seed',1493);
U3 = umod.U;
ok = ok && all(U1 == U3,'all');

% gdata_time, ldata_time
umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata_time' 'Gamma' 'ldata_time'},'src/spin14SIR.c');
umod.data_time = [0 10 20 30 50];
umod.ldata_time = permute([0.2 0.2 0.2 0.2 0.2],[1 3 2]);
umod.gdata_time = [1.5 1.5 1.5 1.5 1.5].*umod.ldata_time(1);
umod.tspan = [0:1:100 1000];
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solve',1,'seed',1493);
U4 = umod.U;
ok = ok && all(U1 == U4,'all');

unix('rm src/spin14SIR.c');

% $$$ figure(1), clf, plot(umod.tspan,umod.U);
% $$$ xline(umod.data_time);
% $$$ legend('S','I','R');

%-------------------------------------------------------------------------------
