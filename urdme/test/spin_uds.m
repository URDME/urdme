function ok = spin_uds(ix)
%SPIN_UDS Tests for UDS solver.

% S. Engblom 2024-05-13 (data_time, ldata_time, gdata_time)
% S. Engblom 2020-02-21 (Revision, Jacobian)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2019-11-17 (Revision, Nreplicas)
% S. Engblom 2017-02-24

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4 @l_spin5};
stests = {'Spin_uds #1' 'Spin_uds #2' 'Spin_uds #3' ...
          'Spin_uds #4' 'Spin_uds #5'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SPIN_UDS (URDME/UDS)', ftests(ix), stests(ix));
unix('rm mexuds*mexrhs.mex*'); % cleanup
unix('rm mexuds*mexjac.mex*');

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
umod.G = sparse(1,4);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin1','solver','uds', ...
             'solverargs',{{'odeopts' odeset('RelTol',1e-6,'AbsTol',1e-4)}});

ok = ok && norm(repmat([1 1 2],1,4)*umod.U-21,inf) < 1e-7;

% Inline propensities.
K = [1 0 0]';
I = [1 2 1]';

vmod = urdme(umod,'propensities','', ...
             'inline_propensities',struct('K',K,'I',I), ...
             'modelname','spin1_inline');
ok = ok && norm(umod.U(:)-vmod.U(:),inf) < 1e-7;

%-------------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Test of reaction-diffusion in four cells with three reactions.
%   X+X->Z, Y+Y->Z and Z->@

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
umod.G = sparse(3,6);

% Output times tspan.
umod.tspan = [0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
umod.sd = ones(1,size(u0,2));

umod = urdme(umod,'propensities','src/spin3','solver','uds');

% Inline propensities.
K = [4 0 0; 4 0 0; 0 1 0]';
I = [1 1 1; 2 2 1; 1 1 3]';
S = sparse(4,3);

vmod = umod;
vmod.propensities = '';
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod.modelname = 'spin3_inline';
vmod = urdme(vmod);
ok = ok && norm(umod.U(:)-vmod.U(:),inf) < 1e-12;

% Nreplicas syntax
vmod.u0 = cat(3,u0+1,u0+2,u0+3);
vmod = urdme(vmod,'modelname','spin3_inline_nreplicas');
for i = 1:size(vmod.u0,3)
  umod.u0 = u0+i;
  umod = urdme(umod,'compile',0,'parse',0);
  ok = ok && norm(umod.U-vmod.U(:,:,i),inf) < 1e-12;
end

%-------------------------------------------------------------------------------
function ok = l_spin3
%L_SPIN3 Five cells with five reactions and four species.
%   D --> A, A+A --> B, B+B --> C, C+C --> A, A+B+C --> D.
%   First reaction turned off for sd == 2.

ok = 1;

% Create the model.
B = [ones(20,2) -3*ones(20,1) ones(20,4)];
B(9:12,1) = 0;
B(5:8,2) = 0;
B(13:16,2) = 0;
B(17:20,3) = 0;
B(9:12,4) = 0;
B(13:16,6) = 0;
D = spdiags(B,[-8 -4 0 4 8 12 16],20,20)';
umod.D = D;

% Initial vector u0.
u0 = [30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0 = u0;

% Stochiometric matrix N in sparse format.
umod.N = sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');

% Dependency graph G in sparse format.
umod.G = sparse([0     0     0     1     1     0     0     0     1; ...
                 1     0     0     0     1     1     0     1     1; ...
                 0     1     0     0     0     1     1     0     1; ...
                 0     0     1     0     0     0     1     1     1; ...
                 1     1     1     0     1     1     1     1     1]);

% Output times tspan.
tspan = 0:10;
umod.tspan = tspan;

% Other input. cell, vol, pos and sd.
umod.vol = ones(1,size(u0,2));
sd = ones(1,size(u0,2));
sd(5) = 2;
umod.sd = sd;

% with/without Jacobian
umod = urdme(umod,'propensities','src/spin6','solver','uds');
U1 = umod.U;
umod = urdme(umod,'propensities','src/spin6','solver','uds', ...
             'solverargs',{{'jacobian' 1}}, ...
             'makeargs',{{'define' '-DJAC_'}});
U2 = umod.U;
u1 = umod.U(:,1);
J1 = mexuds_spin6_mexjac(umod.mexhash,0,u1, ...
                         size(umod.N,2),umod.G, ...
                         umod.vol,[],[],[],[],umod.sd, ...
                         umod.inline_propensities.K, ...
                         umod.inline_propensities.I, ...
                         umod.inline_propensities.S);

% Inline propensities version.
K = [0 1 0;
     2 0 0;
     2 0 0;
     2 0 0]';
I = [1 1 4;
     1 1 1;
     2 2 1;
     3 3 1]';
% D --> A, A+A --> B, B+B --> C, C+C --> A, A+B+C --> D.
S = sparse(zeros(1,4));
S(1,1) = 2;

% using that 4 reactions are inline, the cubic one is first in src/spin6inline.c
vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
% with/without Jacobian
vmod = urdme(vmod,'propensities','src/spin6inline','solver','uds');
V1 = vmod.U;
vmod = urdme(vmod,'propensities','src/spin6inline','solver','uds', ...
             'solverargs',{{'jacobian' 1}}, ...
             'makeargs',{{'define' '-DJAC_'}});
V2 = vmod.U;
J2 = mexuds_spin6inline_mexjac(vmod.mexhash,0,u1, ...
                               size(vmod.N,2),vmod.G, ...
                               vmod.vol,[],[],[],[],vmod.sd, ...
                               vmod.inline_propensities.K, ...
                               vmod.inline_propensities.I, ...
                               vmod.inline_propensities.S);

ok = ok && all([norm(U2-U1,inf) norm(V1-U1,inf) norm(V2-U1,inf)] < 1e-3);
ok = ok && norm(J1-J2,inf) < 1e-12;

% test of rparse/rparse_inline/Jacobian
clear umod;

% fully compiled version incl Jacobian
umod = rparse([],{'D > (sd != 2)*D > A' ...
                  'A+A > A*(A-1)/vol > B' ...
                  'B+B > B*(B-1)/vol > C' ...
                  'C+C > C*(C-1)/vol > A' ...
                  'A+B+C > A*B*C/vol/vol > D'}, ...
              {'A' 'B' 'C' 'D'},{}, ...
              'src/spin6_rparse');
umod = rparse(umod,'jac');

umod.D = D;
umod.u0 = u0;
umod.vol = ones(1,size(u0,2));
umod.sd = sd;
umod.tspan = tspan;
umod = urdme(umod,'solver','uds', ...
             'solverargs',{{'jacobian' 1}}, ...
             'makeargs',{{'define' '-DJAC_'}});
U1_ = umod.U;
J3 = mexuds_spin6_rparse_mexjac(umod.mexhash,0,u1, ...
                                size(umod.N,2),umod.G, ...
                                umod.vol,[],[],[],[],umod.sd, ...
                                umod.inline_propensities.K, ...
                                umod.inline_propensities.I, ...
                                umod.inline_propensities.S);
ok = ok && norm(U2-U1_,inf) < 1e-12; % should be identical this time
ok = ok && norm(J1-J3,inf) < 1e-12;

% mixed version next
clear umod;
umod = rparse([],{'A+B+C > A*B*C/vol/vol > D'}, ...
              {'A' 'B' 'C' 'D'},{},'src/spin6inline_rparse');
umod = rparse_inline(umod,{'D > one > A' 'A+A > two > B' ...
                    'B+B > two > C' 'C+C > two > A'}, ...
                     {'A' 'B' 'C' 'D'},{'one' 1.0 'two' 2.0});
umod.inline_propensities.S(1,1) = 2;
umod = rparse(umod,'jac');
umod.D = D;
umod.u0 = u0;
umod.vol = ones(1,size(u0,2));
umod.sd = sd;
umod.tspan = tspan;
umod = urdme(umod,'solver','uds', ...
             'solverargs',{{'jacobian' 1}}, ...
             'makeargs',{{'define' '-DJAC_'}});
U2_ = umod.U;
J4 = mexuds_spin6inline_rparse_mexjac(umod.mexhash,0,u1, ...
                                      size(umod.N,2),umod.G, ...
                                      umod.vol,[],[],[],[],umod.sd, ...
                                      umod.inline_propensities.K, ...
                                      umod.inline_propensities.I, ...
                                      umod.inline_propensities.S);
ok = ok && norm(V2-U2_,inf) < 1e-12; % should be identical this time
ok = ok && norm(J2-J4,inf) < 1e-12;

unix('rm src/spin6inline_rparse.c'); % cleanup
unix('rm src/spin6_rparse.c');

%-------------------------------------------------------------------------------
function ok = l_spin4
%L_SPIN4 Test of inline propensities. Test of Jacobian.

ok = 1;

k1 = 10; k2 = 12;
nu = 0.001;
mu1 = 0.01; mu2 = 0.09;
ups1 = 0.5; ups2 = 0.7;
umod = rparse_inline([], ...
                     {'@ > k1 > A' '@ > k2 > B' ...
                    'A > mu1 > @' 'B > mu2 > @' ...
                    'A+B > nu > @' ...
                    'B+B > ups1 > A+A' 'A+A > ups2 > B+B'}, ...
                     {'A' 'B'},{'k1' k1 'k2' k2 'nu' nu ...
                    'mu1' mu1 'mu2' mu2 'ups1' ups1 'ups2' ups2});
Mspecies = size(umod.N,1);
Mreactions = size(umod.N,2);

% discretization
Ncells = 20;
x = linspace(0,1,Ncells);
dx = x(2)-x(1);
umod.vol = dx^3*ones(Ncells,1);

% 1D diffusion operator
D = 0.01;
e = ones(Ncells,1)/dx^2;
D = spdiags([e -2*e e],-1:1,Ncells,Ncells);

% reflecting boundary
D(1,1) = -e(1);
D(end,end) = -e(end);

% URDME ordering of species
umod.D = kron(D,eye(2));

% subdomains
umod.sd = ones(1,Ncells);

% initial conditions
umod.u0 = zeros(2,Ncells);

% compile & simulate
umod = urdme(umod,'tspan',0:10,'seed',20200221,'solver','uds', ...
             'solverargs',{{'jacobian',1}});

% exact Jacobian
vol = umod.vol(1);
D = @(A,B)[0 0 mu1 0 nu*B/vol 0 ups2*(2*A-1)/2/vol; ...
           0 0 0 mu2 nu*A/vol ups1*(2*B-1)/2/vol 0]';
u0 = reshape(1:Mspecies*Ncells,Mspecies,Ncells);
D1 = D(u0(1,1),u0(2,1));
[i1,j1,s1] = find(D1);
ii = i1; jj = j1; ss = s1;
for i = 2:Ncells
  Di = D(u0(1,i),u0(2,i));
  i1 = i1+Mreactions;
  j1 = j1+Mspecies;
  [~,~,s1] = find(Di);
  ii = [ii i1]; jj = [jj j1]; ss = [ss s1];
end
Jexact = sparse(ii,jj,ss,Mreactions*Ncells,Mspecies*Ncells);

J = feval([umod.mexname '_mexjac'], ...
          umod.mexhash,0,u0,size(umod.N,2),umod.G, ...
          umod.vol,umod.ldata,umod.gdata, ...
          umod.ldata_time,umod.gdata_time, ...
          umod.sd, ...
          umod.inline_propensities.K, ...
          umod.inline_propensities.I, ...
          umod.inline_propensities.S);
ok = ok && norm(J-Jexact,'fro') < 1e-12;

%-------------------------------------------------------------------------------
function ok = l_spin5
%L_SPIN5 Test of gdata_time, ldata_time.
%   Simple test, copied from spin_ssa(4).

ok = 1;

umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata_time' 'Gamma' 'gdata'},'src/spin5SIR.c');
umod.gdata = 0.2;
umod.data_time = [0 10 20 30 50];
umod.gdata_time = [1.5 3 0 1 5]*umod.gdata;

umod.tspan = 0:1:100;
umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solver','uds');

% switch to scalar ldata_time/ldata-fields instead
vmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata'},'src/spin5SIR.c');
vmod.ldata = 0.2;
vmod.ldata_time = permute([1.5 3 0 1 5]*umod.gdata,[1 3 2]);
vmod.data_time = [0 10 20 30 50];
vmod.tspan = 0:1:100;
vmod.D = sparse(3,3);
vmod.vol = 1e4;
vmod.sd = 1;
vmod.u0 = [vmod.vol-100 100 0]';

vmod = urdme(vmod,'solver','uds');
ok = ok && norm((vmod.U-umod.U)./(1+abs(umod.U)),inf) < 1e-4;

% only ldata_time fields, but in 3 cells
wmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata_time'},'src/spin5SIR.c');
wmod.ldata = [0.2 0.2 0.2];
wmod.ldata_time = permute(repmat([[1.5 3   0   1   5]*umod.gdata; ...
                    0.2 0.2 0.2 0.2 0.2],[1 1 3]),[1 3 2]);
wmod.data_time = [0 10 20 30 50];
wmod.tspan = 0:1:100;
wmod.D = sparse(3,3);
wmod.vol = 1e4;
wmod.sd = 1;
wmod.u0 = [wmod.vol(1)-100 100 0]';
wmod.D = sparse(3*3,3*3);
wmod.vol = [1e4 1e4 1e4];
wmod.sd = [1 1 1];
wmod.u0 = repmat([umod.vol-100 100 0]',1,3);

wmod = urdme(wmod,'solver','uds');
ok = ok && norm((wmod.U(1:3,:)-umod.U)./(1+abs(umod.U)),inf) < 1e-4;

unix('rm src/spin5SIR.c');

% $$$ figure(1), clf, plot(umod.tspan,umod.U);
% $$$ xline(umod.data_time);
% $$$ legend('S','I','R');

%-------------------------------------------------------------------------------
