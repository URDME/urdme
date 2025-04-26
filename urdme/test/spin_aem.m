function ok = spin_aem(ix)
%SPIN_AEM Tests for AEM solver.

% S. Engblom 2024-05-10 (data_time, ldata_time, gdata_time)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2017-02-15 (Minor revision)
% P. Bauer 2013-07-02
% J. Cullhed 2008-06-16

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_spin1 @l_spin2 @l_spin3 @l_spin4 @l_spin5 @l_spin6 @l_spin7 ...
          @l_spin8 @l_spin9 @l_spin10 @l_spin11 @l_spin12 @l_spin13 ...
          @l_spin14 @l_spin15 @l_spin16 @l_spin17};
stests = {'Spin_aem #1' 'Spin_aem #2' 'Spin_aem #3' 'Spin_aem #4' ...
          'Spin_aem #5' 'Spin_aem #6' 'Spin_aem #7' 'Spin_aem #8' ...
          'Spin_aem #9' 'Spin_aem #10' 'Spin_aem #11' 'Spin_aem #12' ...
          'Spin_aem #13' 'Spin_aem #14' 'Spin_aem #15' ...
          'Spin_aem #16' 'Spin_aem #17'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SPIN_AEM (URDME/AEM)', ftests(ix), stests(ix));
unix('rm mexaem*.mex*'); % cleanup

%-------------------------------------------------------------------------------
function ok=l_spin1
%L_SPIN1 Test of reaction-diffusion in four cells with one reaction.
%   X+Y --> Z

ok=1;

% Create a simple model.
B=[ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2)=0;
B(7:9,4)=0;

D=spdiags(B,[-6 -3 0 3 6],12,12);
umod.D=D';

% Initial vector u0.
u0=[3 3 0; 3 3 0; 3 3 0; 2 1 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1 -1 1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1 1 1]);

% Output times tspan.
umod.tspan=[0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin1','report',0,'solver','aem','seed',123);

xx=umod.U(:,length(umod.tspan));
ok=ok&&sum([xx(1) xx(4) xx(7) xx(10)])==1;
ok=ok&&sum([xx(2) xx(5) xx(8) xx(11)])==0;
ok=ok&&sum([xx(3) xx(6) xx(9) xx(12)])==10;
ok=ok&&all(umod.U(:)>=0);

% Inline propensities.
K=[1 0 0]';
I=[1 2 1]';

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod.compile = 0; % don't re-compile
vmod = urdme(vmod);
ok = ok && all(umod.U(:) == vmod.U(:));

% 1st test of the AEM-property
wmod = vmod;
K1 = 1.0005*K; % *** unclear: apparently this perturbation can be
               % larger on some platforms...?
wmod.inline_propensities = struct('K',K1,'I',I);
wmod = urdme(wmod);
ok = ok && all(vmod.U(:) == wmod.U(:));

% Nreplicas with different seeds
seed = [123 456 765 433 223];
wmod.inline_propensities = struct('K',K,'I',I);
wmod.u0 = repmat(u0,1,1,numel(seed));
wmod.seed = seed;
wmod = urdme(wmod);
for i = 1:numel(seed)
  vmod.seed = seed(i);
  vmod = urdme(vmod);
  ok = ok && norm(vmod.U-wmod.U(:,:,i),'fro') == 0;
end

% very small perturbation
seed = [123 456 765 433 223];
wmod.inline_propensities = struct('K',1.0001*K,'I',I);
wmod.u0 = repmat(u0,1,1,numel(seed));
wmod.seed = seed;
wmod = urdme(wmod);
for i = 1:numel(seed)
  vmod.seed = seed(i);
  vmod = urdme(vmod);
  ok = ok && norm(vmod.U-wmod.U(:,:,i),'fro') == 0;
end

% a little larger
seed = [123 456 433 223];
wmod.inline_propensities = struct('K',1.0005*K,'I',I);
wmod.u0 = repmat(u0,1,1,numel(seed));
wmod.seed = seed;
wmod = urdme(wmod);
for i = 1:numel(seed)
  vmod.seed = seed(i);
  vmod = urdme(vmod);
  ok = ok && norm(vmod.U-wmod.U(:,:,i),'fro') == 0;
end

% a little larger
seed = [456];
wmod.inline_propensities = struct('K',1.001*K,'I',I);
wmod.u0 = repmat(u0,1,1,numel(seed));
wmod.seed = seed;
wmod = urdme(wmod);
for i = 1:numel(seed)
  vmod.seed = seed(i);
  vmod = urdme(vmod);
  ok = ok && norm(vmod.U-wmod.U(:,:,i),'fro') == 0;
end

%-------------------------------------------------------------------------------
function ok=l_spin2
%L_SPIN2 Test of reaction-diffusion in four cells with two reactions.
% X+X->Z and Y+Y->Z
ok=1;

% Create a simple model.
B=[ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2)=0;
B(7:9,4)=0;
D=spdiags(B,[-6 -3 0 3 6],12,12);
umod.D=D';

% Initial vector u0.
u0=[10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-2 -0 1; -0 -2 1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 0 0 1 0; 0 1 0 0 1]);

% Output times tspan.
umod.tspan=[0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin2','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&sum([xx(1) xx(4) xx(7) xx(10)])==0;
ok=ok&&sum([xx(2) xx(5) xx(8) xx(11)])==0;
ok=ok&&sum([xx(3) xx(6) xx(9) xx(12)])==50;
ok=ok&&all(umod.U(:)>=0);

%-------------------------------------------------------------------------------
function ok=l_spin3
%L_SPIN3 Test of reaction-diffusion in four cells with three reactions.
% X+X->Z, Y+Y->Z and Y->Z
ok=1;
 
% Create a simple model. 
B=[ones(12,2) -2*ones(12,1) ones(12,2)]; 
B(4:6,2)=0; 
B(7:9,4)=0;
D=spdiags(B,[-6 -3 0 3 6],12,12); 
umod.D=D';

% Initial vector u0. 
u0=[10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0=u0;
 
% Stochiometric matrix N in sparse format.
umod.N=sparse([-2 -0 1; 0 -2 1; 0 0 -1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 0 0 1 0 0; 0 1 0 0 1 0; 0 0 1 1 1 1]);

% Output times tspan.
umod.tspan=[0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin3','report',0,'solver','aem','seed',29);

xx=umod.U(:,length(umod.tspan));
ok=ok&&sum([xx(1) xx(4) xx(7) xx(10)])==0;
ok=ok&&sum([xx(2) xx(5) xx(8) xx(11)])==0;
ok=ok&&sum([xx(3) xx(6) xx(9) xx(12)])==0;
ok=ok&&all(umod.U(:)>=0);

% Inline propensities.
K=[4 0 0; 4 0 0; 0 1 0]';
I=[1 1 1; 2 2 1; 1 1 3]';
S=sparse(4,3);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok=l_spin4
%L_SPIN4 Test of reaction-diffusion in four cells with three reactions.
% X+X->Z and Y+Y->Z
ok=1;

% Create a simple model. 
B=[ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2)=0;
B(7:9,4)=0;
D=spdiags(B,[-6 -3 0 3 6],12,12);
D(3:3:end,:)=0;
D(:,3:3:end)=0;
umod.D=D';

% Initial vector u0. 
u0=[10 37 0; 10 3 0; 10 10 0; 20 0 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-2 0 1; 0 -2 1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 0 0 1 0; 0 1 0 0 1]);

% Output times tspan.
umod.tspan=[0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
sd=ones(1,size(u0,2));
sd(1)=2;
umod.sd=sd;

umod=urdme(umod,'propensities','src/spin4','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&xx(3)==50;
ok=ok&&all(umod.U(:)>=0);

%-------------------------------------------------------------------------------
function ok=l_spin5
%L_SPIN5 Test of reaction-diffusion in four cells with one reaction and
%diffusion in all cells except cell one.
% X+Y->Z
ok=1;

% Create a simple model.
B=[ones(12,2) -2*ones(12,1) ones(12,2)];
B(4:6,2)=0;
B(1:3,3)=0;
B(4:9,4)=0;
B(7:9,5)=0;
D=spdiags(B,[-6 -3 0 3 6],12,12);
umod.D=D';

% Initial vector u0.
u0=[3 3 0; 3 3 0; 3 3 0; 2 1 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1 -1 1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1 0 1]);

% Output times tspan.
umod.tspan=[0:10 70:10:200];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));
umod=urdme(umod,'propensities','src/spin5','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&xx(1)==1&&xx(3)==10;
ok=ok&&all(umod.U(:)>=0);

%-------------------------------------------------------------------------------
function ok=l_spin6
%L_SPIN6 Test of reaction-diffusion in five cells with five reactions and four
%species.
% D -> A, A+A -> B, B+B -> C, C+C -> A, A+B+C -> D.
ok=1;

% Create the model.
B=[ones(20,2) -3*ones(20,1) ones(20,4)];
B(9:12,1)=0;
B(5:8,2)=0;
B(13:16,2)=0;
B(17:20,3)=0;
B(9:12,4)=0;
B(13:16,6)=0;
D=spdiags(B,[-8 -4 0 4 8 12 16],20,20);
umod.D=D';

% Initial vector u0.
u0=[30 30 0 3; 32 31 3 9; 3 3 10 40; 2 1 0 7; 10 10 10 10]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([1 0 0 -1; -2 1 0 0; 0 -2 1 0; 1 0 -2 0; -1 -1 -1 1]');

% Dependency graph G in sparse format.
umod.G=sparse([0 0 0 1 1 0 0 0 1;
               1 0 0 0 1 1 0 1 1;
               0 1 0 0 0 1 1 0 1;
               0 0 1 0 0 0 1 1 1;
               1 1 1 0 1 1 1 1 1]);

% Output times tspan.
umod.tspan=[0:1:100];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
sd=ones(1,size(u0,2));
sd(5)=2;
umod.sd=sd;

umod=urdme(umod,'propensities','src/spin6','seed',20120303, ...
           'report',0,'solver','aem');
% *** platform dependent:
%ok=ok&&sum(umod.U(:,end))==72; % or 70? 
ok=ok&&all(umod.U(:)>=0);

% Inline propensities.
K=[0 1 0;
   2 0 0;
   2 0 0;
   2 0 0]';

I=[1 1 4;
   1 1 1;
   2 2 1;
   3 3 1]';
% D -> A, A+A -> B, B+B -> C, C+C -> A, A+B+C -> D.
S=sparse(zeros(1,4));
S(1,1)=2;

% using that 4 reactions are inline, the cubic one is first in src/spin6inline.c
vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'propensities','src/spin6inline', ...
             'seed',20120303);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok=l_spin7
%L_SPIN7 Test of reaction-diffusion in six cells with one reaction and three
%species.
% X+Y->Z, X->Z, Y->Z
ok=1;

% Create the model.
B=[-1*ones(18,1) ones(18,1)];
B(16:18,1)=0;
D=spdiags(B,[0 3],18,18);
umod.D=D';

% Initial vector u0.
u0=[1000 1000 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1 -1 1; -1 1 0; 1 -1 0]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1 0 1 1 1;
          1 0 0 1 1 1;
          0 1 0 1 1 1]);

% Output times tspan.
umod.tspan=[0:0.3:101];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin7','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&xx(18)==1000;
ok=ok&&all(umod.U(:)>=0);

%-------------------------------------------------------------------------------
function ok=l_spin8
%L_SPIN8 Test of reaction-diffusion in four cells with three reactions and three
%species. Diffusion is in opposit directions for x,w and y,z.
% W+X->Z, 100*(W+Y)->Z, X+Y->Z
ok=1;

% Create the model.
B1=[0 1 1 0 0 1 1 0 0 1 1 0 0 0 0 0]';
B2=[-1 0 0 -1*ones(1,9) 0 -1 -1 0]';
B3=[0 0 0 0 1 0 0 1 1 0 0 1 1 0 0 1]';

B=[B1 B2 B3];
D=spdiags(B,[-4 0 4],16,16);
umod.D=D';

% Initial vector u0.
u0=[100 0 0 100; 0 0 0 0; 0 0 0 0; 0 1000 100 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1 -1 0 1; 1 -1 -1 0; -1 1 -1 0]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1 0 0 1 1 1;
          0 1 1 0 1 1 1;
          1 0 1 0 1 1 1]);

% Output times tspan.
umod.tspan=0:1:100;

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin8','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&xx(2)==700&&xx(end)==300;
ok=ok&&all(umod.U(:)>=0);
%-------------------------------------------------------------------------------
function ok=l_spin9
%L_SPIN9 Test of reaction-diffusion in four cells with one reaction and one
%species. This is a test of the trivial case.
% Y -> Z
ok=1;

% Create the model.
D=sparse([-2 1 1 0; 1 -2 0 1; 1 0 -2 1; 0 1 1 -2]);
umod.D=D';

% Initial vector u0.
u0=[1; 100; 3; 5]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1]);

% Output times tspan.
umod.tspan=[0:0.321:101];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod=urdme(umod,'propensities','src/spin9','report',0,'solver','aem');

xx=umod.U(:,length(umod.tspan));
ok=ok&&any(xx==0);
ok=ok&&all(umod.U(:)>=0);
%-------------------------------------------------------------------------------
function ok=l_spin10
%L_SPIN10 Test of random seed. Compare two models with the same random seed
%and see if the result is the same.
% X+Y->Z, X->Y and Y->Z
ok=1;

% Create the model.
B=[-1*ones(18,1) ones(18,1)];
B(16:18,1)=0;
D=spdiags(B,[0 3],18,18);
umod.D=D';

% Initial vector u0.
u0=[1000 1000 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0]';
umod.u0=u0;

% Stochiometric matrix N in sparse format.
umod.N=sparse([-1 -1 1; -1 1 0; 1 -1 0]');

% Dependency graph G in sparse format.
umod.G=sparse([1 1 0 1 1 1;
          1 0 0 1 1 1;
          0 1 0 1 1 1]);

% Output times tspan.
umod.tspan=[0:0.3:101];

% Other input. cell, vol, pos and sd.
umod.vol=ones(1,size(u0,2));
umod.sd=ones(1,size(u0,2));

umod1=urdme(umod,'propensities','src/spin10','seed',20120329,'report',0, ...
            'solver','aem');
umod2=urdme(umod,'propensities','src/spin10','seed',20120329,'report',0, ...
            'solver','aem');

ok=ok&&all(all(umod1.U==umod2.U));
ok=ok&&all(umod1.U(:)>=0);

%-------------------------------------------------------------------------------
function ok=l_spin11
%L_SPIN11: Testing of correct reaction and diffusion waiting times;
%   1) 2 Reaction channels
%   2) 2 Diffusion channels
%   Channel rates = {1,10}

ok = 1;

% 1) Diffusion case

% Create the model.
D = sparse([-1  1;
           10 -10]);
umod1.D = D';

% Initial vector u0.
u0 = [10; 0]';
umod1.u0 = u0;

% Stochiometric matrix N in sparse format.
umod1.N = sparse(zeros(1,0));

% Dependency graph G in sparse format.
umod1.G = sparse(zeros(0,1));

% Output times tspan.
umod1.tspan = 0:.05:1;

% Other input. cell, vol, pos and sd.
umod1.vol = ones(1,size(u0,2));
umod1.sd = ones(1,size(u0,2));

umod1 = urdme(umod1,'propensities','src/spin11','report',0, ...
              'solver','aem','seed',1234);

% 1) Reaction case

% Geometry - single cell
D = sparse([0 0 ; 0 0]);
umod2.D = D';

% Initial vector u0.
u0 = [10 0]';
umod2.u0 = u0;

% Stochiometric matrix N in sparse format.
umod2.N = sparse([-1 1; 1 -1]');

% Dependency graph G in sparse format.
umod2.G = sparse([1 1 1 1; 1 1 1 1]);

% Output times tspan.
umod2.tspan = umod1.tspan;

% Other input. cell, vol, pos and sd.
umod2.vol = ones(1,size(u0,2));
umod2.sd = ones(1,size(u0,2));

umod2 = urdme(umod2,'propensities','src/spin11','report',0, ...
              'solver','aem','seed',1234);

ok = ok && isequal(umod1.U,umod2.U);

%-------------------------------------------------------------------------------
function ok=l_spin12
%L_SPIN11: Testing of correct reaction and diffusion waiting times;
%   1) 3 Reaction channels
%   2) 3 Diffusion channels
%   Channel rates = {2,5,8}

ok = 1;

% 1) Diffusion case

% Create the model.
D = sparse([-2    2    0;
            0   -5    5;
            8    0   -8]);
umod1.D = D';

% Initial vector u0.
u0 = [10 10 10];
umod1.u0 = u0;

% Stochiometric matrix N in sparse format.
umod1.N = sparse(zeros(1,0));

% Dependency graph G in sparse format.
umod1.G = sparse(zeros(0,1));

% Output times tspan.
umod1.tspan = 0:1:100;

% Other input. cell, vol, pos and sd.
umod1.vol = ones(1,size(u0,2));
umod1.sd = ones(1,size(u0,2));

umod1 = urdme(umod1,'propensities','src/spin12','report',0, ...
              'solver','aem','seed',1234567);

% 2) Reaction case

% Geometry - single cell
D = sparse([0 0 0;
            0 0 0;
            0 0 0]);
umod2.D = D';

% Initial vector u0.
u0 = [10 10 10]';
umod2.u0 = u0;

% Stochiometric matrix N in sparse format.
umod2.N = sparse([-1  1  0;
                  0 -1  1;
                  1  0 -1]');

% Dependency graph G in sparse format.
umod2.G = sparse([0 0 0 1 0 1; 
                  0 0 0 1 1 0;
                  0 0 0 0 1 1]);
                     
% Output times tspan.
umod2.tspan = umod1.tspan;

% Other input. cell, vol, pos and sd.
umod2.vol = ones(1,size(u0,2));
umod2.sd = ones(1,size(u0,2));

umod2 = urdme(umod2,'propensities','src/spin12','report',0, ...
              'solver','aem','seed',1234567);

ok = ok && isequal(umod1.U,umod2.U);

%-------------------------------------------------------------------------------
function ok=l_spin13
%L_SPIN13: Test of correct reaction and diffusion waiting on 4 domains
%   1) 4 Reaction channels
%   2) 4 Diffusion channels
%   Channel rates = {2,5,8,17}
%   Reaction seeds are given away before diffusion seeds

ok = 1;

% 1) 4 Reaction channels
% A --> B, B --> C, C --> D, D --> A

% Create the model.
D = sparse([0 0 0 0; 
            0 0 0 0; 
            0 0 0 0;
            0 0 0 0]);
umod1.D = D';

% Initial vector u0.
u0 = [100 100 100 100]';
umod1.u0 = u0;

% Stochiometric matrix N in sparse format.
umod1.N = sparse([-1  1  0  0; 
                  0 -1  1  0;
                  0  0 -1  1;
                  1  0  0 -1]');

% Dependency graph G in sparse format.
umod1.G = sparse([0 0 0 0 1 0 0 1; 
                  0 0 0 0 1 1 0 0;
                  0 0 0 0 0 1 1 0;
                  0 0 0 0 0 0 1 1]);

% Output times tspan.
umod1.tspan = 0:1:100;

% Other input. cell, vol, pos and sd.
umod1.vol = ones(1,size(u0,2));
umod1.sd = ones(1,size(u0,2));

umod1 = urdme(umod1,'propensities','src/spin13','report',0, ...
              'solver','aem','seed',1234);

% 2) 4 Diffusion channels
% C1 -> C2 -> C3 -> C4 -> C1

D = sparse([-2    2    0   0;
            0   -5    5   0;
            0    0   -8   8;
            17    0    0 -17]);
umod2.D = D';

% Initial vector u0.
u0 = [100 100 100 100];
umod2.u0 = u0;

% Stochiometric matrix N in sparse format.
umod2.N = sparse(zeros(1,0));

% Dependency graph G in sparse format.
umod2.G = sparse(zeros(0,1));

% Output times tspan.
umod2.tspan = umod1.tspan;

% Other input. cell, vol, pos and sd.
umod2.vol = ones(1,size(u0,2));
umod2.sd = ones(1,size(u0,2));

umod2 = urdme(umod2,'propensities','src/spin13','report',0, ...
              'solver','aem','seed',1234);

ok = ok && isequal(umod1.U,umod2.U);

%-------------------------------------------------------------------------------
function ok=l_spin14
%L_SPIN14: Test of correct reaction and diffusion waiting on 3 domains
%  1) 4 Reaction channels
%  2) 4 Diffusion channels
%  Channel rates = {2,5,8,17}
%  Reaction seeds are given away before diffusion seeds

ok = 1;

% 1) 4 Reaction channels, 0d
% A --> B, B --> A, B --> C, C --> B

% Create the model.
D = sparse([0 0 0; 0 0 0; 0 0 0]);
umod1.D = D';

% Initial vector u0.
u0 = [10 20 5];
umod1.u0 = u0';

% Stochiometric matrix N in sparse format.
umod1.N = sparse([-1  1  0; 
                  1 -1  0;
                  0 -1  1;
                  0  1 -1]');

% Dependency graph G in sparse format.
umod1.G = sparse([0 0 0 1 1 0 0; 
                  0 0 0 1 1 1 1;
                  0 0 0 1 1 1 1;
                  0 0 0 0 0 1 1]);

% Output times tspan.
umod1.tspan = 0:.01:.1;

% Other input. cell, vol, pos and sd.
umod1.vol = 1;
umod1.sd  = 1;
umod1.ldata = 1;

umod1 = urdme(umod1,'propensities','src/spin14','report',0, ...
              'solver','aem','seed',1234);

% 2) 4 Diffusion channels
D = sparse([-2    2    0; ...
            5  -13    8; ...
            0   17   -17]);
umod2.D = D';

% Initial vector u0.
umod2.u0 = u0;

% Stochiometric matrix N in sparse format.
umod2.N = sparse(zeros(1,0));

% Dependency graph G in sparse format.
umod2.G = sparse(zeros(0,1));

% Output times tspan.
umod2.tspan = umod1.tspan;

% Other input. cell, vol, pos and sd.
umod2.vol = ones(1,size(u0,2));
umod2.sd = ones(1,size(u0,2));

umod2 = urdme(umod2,'propensities','src/spin14','report',0, ...
              'solver','aem','seed',1234);

ok = ok && isequal(umod1.U,umod2.U);

%-------------------------------------------------------------------------------
function ok = l_spin15
%L_SPIN15 Linear birth-death model, parameters passed in gdata.
%   Test of the AEM-property.

ok = 1;

% diffusion
D = sparse([-1 1 0; 1 -2 1; 0 1 -1]);
umod.D = D';

% initial vector u0
umod.u0 = zeros(1,3);

% stochiometric matrix
umod.N = sparse([1 -1]); 

% dependency graph
umod.G = sparse([0 0 0; ...
                 1 1 1]);

% Output times tspan.
umod.tspan = linspace(0,1,100);

% vol, sd, data
umod.vol = ones(1,3);
umod.sd = ones(1,3);
umod.ldata = zeros(0,3);
umod.gdata = [1e3 10]; % gdata = [birth death], koefficients of propensities
umod.data_time = 0;
umod.ldata_time = zeros(0,0,3);
umod.gdata_time = zeros(0,0);
umod.propensities = 'src/spin15.c';

% inline propensities
K = [0 0; ...
     0 umod.gdata(2); ...
     umod.gdata(1) 0];
I = ones(3,2);
umod.makeargs = {};
umod.solverargs = {};

umod.solver = 'aem';
umod.seed = 1729;
umod.inline_propensities = struct('K',[],'I',[],'S',[]);
umod.report = 0;
umod.solve = 1;
umod.parse = 0;
umod.modelname = 'test_spin15';
umod.mexhash = 0xFF11FF00;
umod.mexname = 'mexaem_test_spin15';
umod.compile = 1;

umod = urdme(umod); % (note: not parsed)
umod.compile = 0;
U1 = umod.U;

umod.gdata(1) = umod.gdata(1)*1.001;
umod = urdme(umod);
U2 = umod.U;
ok = ok && norm(U1-U2,'fro')/norm(U1,'fro') < 0.02;

umod.gdata(2) = umod.gdata(2)*1.001;
umod = urdme(umod);
U3 = umod.U;
ok = ok && norm(U1-U3,'fro')/norm(U1,'fro') < 0.02;

umod.inline_propensities = struct('K',K,'I',I,'S',[]);
umod = urdme(umod);
V1 = umod.U;
ok = ok && norm(U1-V1,'fro')/norm(U1,'fro') < 0.02;

umod.inline_propensities.K(:,1) = umod.inline_propensities.K(:,1)*1.001;
umod = urdme(umod);
V2 = umod.U;
ok = ok && norm(V1-V2,'fro')/norm(V1,'fro') < 0.02;

umod.inline_propensities.K(:,2) = umod.inline_propensities.K(:,2)*1.001;
umod = urdme(umod);
V3 = umod.U;
ok = ok && norm(V1-V3,'fro')/norm(V1,'fro') < 0.02;

%-------------------------------------------------------------------------------
function ok = l_spin16
%L_SPIN16 Test of gdata_time, ldata_time.
%   Simple test, copied from spin_ssa(4).

ok = 1;

% gdata_time/gdata-fields
umod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'gdata_time' 'Gamma' 'gdata'},'src/spin16SIR.c');
umod.gdata = 0.2;
umod.data_time = [0 10 20 30 50];
umod.gdata_time = [1.5 3 0 1 5]*umod.gdata;
umod.tspan = 0:1:100;

umod.D = sparse(3,3);
umod.vol = 1e4;
umod.sd = 1;
umod.u0 = [umod.vol-100 100 0]';
umod = urdme(umod,'solver','aem','seed',19443,'solve',1);

% switch to scalar ldata_time/ldata-fields instead
vmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata'},'src/spin16SIR.c');
vmod.ldata = 0.2;
vmod.data_time = [0 10 20 30 50];
vmod.ldata_time = permute([1.5 3 0 1 5]*umod.gdata,[1 3 2]);

vmod.tspan = 0:1:100;
vmod.D = sparse(3,3);
vmod.vol = 1e4;
vmod.sd = 1;
vmod.u0 = [vmod.vol-100 100 0]';
vmod = urdme(vmod,'solver','aem','seed',19443);

ok = ok && all(umod.U == vmod.U,'all');

% only ldata_time fields, but in 3 cells
wmod = rparse([],{'S+I > Beta*S*I/vol > I+I' ...
                  'I > Gamma*I > R'},{'S' 'I' 'R'}, ...
              {'Beta' 'ldata_time' 'Gamma' 'ldata_time'},'src/spin16SIR.c');
wmod.ldata = [0.2 0.2 0.2];
wmod.data_time = [0 10 20 30 50];
wmod.ldata_time = permute(repmat([[1.5 3   0   1   5]*umod.gdata; ...
                    0.2 0.2 0.2 0.2 0.2],[1 1 3]),[1 3 2]);
wmod.tspan = 0:1:100;
wmod.D = sparse(3*3,3*3);
wmod.vol = [1e4 1e4 1e4];
wmod.sd = [1 1 1];
wmod.u0 = repmat([wmod.vol(1)-100 100 0]',1,3);

wmod = urdme(wmod,'solver','aem','seed',19443);
ok = ok && all(wmod.U(1:3,:) == umod.U,'all');

unix('rm src/spin16SIR.c');

% $$$ figure(1), clf, plot(umod.tspan,umod.U);
% $$$ xline(umod.data_time);
% $$$ legend('S','I','R');

%-------------------------------------------------------------------------------
function ok = l_spin17
%L_SPIN17 Rectangular domain with diffusion into another species.
%   Simple test, copied from scrub(5).

ok = 1;

% fetch discretization
N = 5;
[P,E,T] = basic_mesh(2,N);

% assemble diffusion operator and voxel volume vector
umod = pde2urdme(P,T,{1});

% reactions
umod.sd(:) = 1;
Ncells = numel(umod.sd);

rates = {'b' 2 'c' 1};
species = {'A' 'B' 'C'};

% transitions and rates
r = cell(1,2);
r{1} = 'B > b > C';
r{2} = 'C > c > @';

% check that the special diffusion is noticed:
warning('off','rparse_inline:ghost_species');
% sort it out
umod = rparse_inline(umod,r,species,rates);
[~,wid] = lastwarn;
warning('on','rparse_inline:ghost_species');
ok = ok && strcmpi(wid,'rparse_inline:ghost_species');
Mspecies = size(umod.N,1);

% define the diffusion matrix
L = umod.D/min(abs(diag(umod.D)));

% the negative (outgoing) part of L concerns species A, the positive
% (incoming) part of L concerns B, i.e., A diffuses into a neighboring
% voxel, transforming into a B at the same time
ixA = strcmp(umod.private.rpinl.Species,'A');
ixB = strcmp(umod.private.rpinl.Species,'B');
% perhaps not 100% obvious, but this is how it is done:
umod.D = kron(L.*(L < 0),diag(sparse(ixA)))+ ...
         kron(L.*(L > 0),sparse(find(ixB),find(ixA),1,Mspecies,Mspecies));

umod.u0 = [100*ones(1,Ncells); zeros(2,Ncells)];
% so A diffuses into a B which transform into a C which is destroyed...
umod.tspan = [0:1:10 100];
umod.seed = 240411;

% check that the special diffusion is noticed:
warning('off','urdme_validate_model:species_diffusion');
urdme_validate_model(umod,'d');
[~,wid] = lastwarn;
warning('on','urdme_validate_model:species_diffusion');
ok = ok && strcmpi(wid,'urdme_validate_model:species_diffusion');

umod = urdme(umod,'solver','aem');
% ...after some time nothing remains:
ok = ok && all(reshape(umod.U(:,end),3,[]) == 0,'all');

%-------------------------------------------------------------------------------
