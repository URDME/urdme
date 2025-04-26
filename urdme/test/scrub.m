function ok = scrub(ix)
%SCRUB Tests for URDME.

% S. Engblom 2024-04-11 (fresh set of scrub*.mat and added advanced diffusion example)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2017-02-18 (Revision, inline propensities)
% S. Engblom 2017-02-15 (Minor revision)
% P. Bauer 2012-08-31
% J. Cullhed 2008-08-08

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_scrub0 @l_scrub1 @l_scrub2 @l_scrub3 @l_scrub4 @l_scrub5};
stests = {'Scrub #0' 'Scrub #1' 'Scrub #2' 'Scrub #3' 'Scrub #4' 'Scrub #5'};

if nargin == 0
  ix = 1:size(ftests,2);
else
  ix = ix+1; % zero-offset, this scrub only
end

ok = runtest('SCRUB (URDME)', ftests(ix), stests(ix));
unix('rm mexnsm*.mex*');
unix('rm mexaem*.mex*');
unix('rm mexuds*.mex*');
unix('rm mexssa*.mex*');

%-------------------------------------------------------------------------------
function ok = l_scrub0
%L_SCRUB0 Run help texts, should not crash...

ok = 1;

% NSM
str = help('nsm');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1)); % eval what is in between

% AEM
str = help('aem');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% SSA, but this actually runs the UDS-solver too:
str = help('ssa');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% UDS, but this actually runs the NSM-solver too;
str = help('uds');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% DLCM:
str  = help('dlcm');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% DLCM2URDME:
str  = help('dlcm2urdme');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% DT_OPERATORS:
str  = help('dt_operators');
i = findstr(str,'Example:');
j = findstr(str,'See also');
eval(str(i+8:j-1));

% RPARSE/RPARSE_INLINE
str = help('rparse');
i = findstr(str,'Examples:');
j = findstr(str,'See also');
eval(str(i+9:j-1));
unix('rm bimol.c');
str = help('rparse_inline');
i = findstr(str,'Examples:');
j = findstr(str,'See also');
eval(str(i+9:j-1));

% SEQEXPAND
str = help('seqexpand');
i = findstr(str,'Examples:');
eval(str(i+9:end));

% URDME_VALIDATE_MODEL
str = help('urdme_validate_model');
i = findstr(str,'Example:');
j = findstr(str,'See also');
try
  eval(str(i+9:j-1));
  ok = 0;
catch
  msg = 'Reaction #2 may produce negative states.';
  le = lasterr;
  ok = ok && ...
       strcmp(le(end-numel(msg)+1:end),msg);
end
unix('rm test.c');

close all;

%-------------------------------------------------------------------------------
function ok = l_scrub1
%L_SCRUB1 Test of reaction-diffusion in 4291 cells with one reaction.

ok = 1; 
load data/scrub1
umod = urdme(umod,'propensities','src/scrub1','report',0,'seed',1234);

ok = ok && sum(umod.U(1:2:end,end)) == 0 && ...
     sum(umod.U(2:2:end,end)) == 1000;
ok = ok && all(umod.U(:) >= 0);

% compare using inline propensities
K = [1.9293e-14 0 0]';
I = [1 2 1]';

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod = urdme(vmod,'seed',1234,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok = l_scrub2
%L_SCRUB2 Test of reaction-diffusion in 4291 cells with two reactions.

ok = 1; 
load data/scrub2
umod = urdme(umod,'propensities','src/scrub2','seed',345345);

ok = ok && all(sum(reshape(umod.U(:,end),4,[]),2)==[0 6000 0 1000]');
ok = ok && all(umod.U(:) >= 0);

% compare using inline propensities
K = [1.9293e-14 0 0; ...
     3.8585e-14 0 0]';
I = [1 2 1; ...
     3 4 1]';
    
vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I);
vmod = urdme(vmod,'seed',345345,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok = l_scrub3
%L_SCRUB3 Test of reaction-diffusion in 700 cells with four reactions.

ok = 1;
load data/scrub3
umod = urdme(umod,'propensities','src/scrub3','seed',678678);

ok = ok && sum(umod.U(4:5:end,end)) == 93;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [0 0 1e4; ...
     0 0 1e4; ...
     0.01 0 0; ...
     0.01 0 0]';
I = [1 1 1; ...
     1 1 1; ...
     1 2 1; ...
     3 5 1]';
S = sparse([1 1 0 0]);

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'seed',678678,'compile',0);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok = l_scrub4
%L_SCRUB4 Test of reaction-diffusion in 1509 cells with three reactions.

ok = 1;
load data/scrub4
umod = urdme(umod,'propensities','src/scrub4','report',0,'seed',20120330);

ok = ok && sum(umod.U(3:3:end,end))==443;
ok = ok && all(umod.U(:) >= 0);

% Inline propensities.
K = [0 0 1e3; ...
     0 1 0; ...
     0 1 0]';

I = [1 1 1; ...
     1 1 1; ...
     1 1 2]';
S = sparse([1 2 4; 1 2 3; 1 3 4]');

vmod = umod;
vmod.inline_propensities = struct('K',K,'I',I,'S',S);
vmod = urdme(vmod,'seed',20120330);
ok = ok && all(umod.U(:) == vmod.U(:));

%-------------------------------------------------------------------------------
function ok = l_scrub5
%L_SCRUB5 Model in a rectangular domain with diffusion into another species.

ok = 1;

% fetch discretization
N = 10;
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
umod.tspan = 0:10:100;
umod.seed = 240411;

% check that the special diffusion is noticed:
warning('off','urdme_validate_model:species_diffusion');
urdme_validate_model(umod,'d');
[~,wid] = lastwarn;
warning('on','urdme_validate_model:species_diffusion');
ok = ok && strcmpi(wid,'urdme_validate_model:species_diffusion');

umod = urdme(umod,'solver','nsm','tspan',0:10:100,'seed',240411);
% ...after some time nothing remains:
ok = ok && all(reshape(umod.U(:,end),3,[]) == 0,'all');

%-------------------------------------------------------------------------------
