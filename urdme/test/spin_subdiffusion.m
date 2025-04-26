function ok = spin_subdiffusion(ix)
%SPIN_SUBDIFFUSION Tests for SUBDIFFUSION workflow.

% S. Engblom 2017-05-16

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

ftests = {@l_spin1 @l_spin2};
stests = {'Spin_subdiffusion #1' 'Spin_subdiffusion #2'};
if nargin == 0, ix = 1:size(ftests,2); end
ok = runtest('SPIN_SUBDIFFUSION (WORKFLOWS/SUBDIFFUSION)', ...
             ftests(ix),stests(ix));

%-------------------------------------------------------------------------------
function ok = l_spin1
%L_SPIN1 Test of SUBDMATRIX.

ok = 1;

% mock internal states transfer matrix
A = reshape(12:15,2,2); % 2 states
A = A-diag(sum(A,1));
gamma = [1 15];
D = full(spdiags(reshape(1:9,3,3),[-1 0 1],3,3)); % 3 voxels
D = D-diag(sum(D,1));
D = kron(D,diag(1:5)); % 5 species

% select all species
select = [1 1 1 1 1];
SD = full(subdmatrix(A,gamma,D,logical(select)));
SD2 = full(subdmatrix(A,gamma,D));

ok = ok && norm(SD-SD2,'fro') < 1e-12;
ok = ok && norm(sum(SD,1),inf) < 1e-12;

% select a few species
select = [1 0 1 0 1];
SD = full(subdmatrix(A,gamma,D,logical(select)));
state1 = cumsum([1 1+select*(size(A,1)-1)]);
Mspecies = state1(end)-1;
[ii,jj] = find(SD);
from_vol = ceil(ii/Mspecies);
from_spec = mod(ii,Mspecies);
to_vol = ceil(jj/Mspecies);
to_spec = mod(jj,Mspecies);
ok = ok && sum(from_vol ~= to_vol & from_spec ~= to_spec) == 0;
ok = ok && norm(sum(SD,1),inf) < 1e-12;

% slightly larger test

A = reshape(1:16,4,4); % 4 states
A = A-diag(sum(A,1));
gamma = [1 15 2 16];
D = full(spdiags(reshape(1:25,5,5),[-1 0 1],5,5)); % 5 voxels
D = D-diag(sum(D,1));
D = kron(D,diag(1:5)); % 5 species

select = [1 0 1 1 0];
SD = full(subdmatrix(A,gamma,D,logical(select)));
state1 = cumsum([1 1+select*(size(A,1)-1)]);
Mspecies = state1(end)-1;
[ii,jj] = find(SD);
from_vol = ceil(ii/Mspecies);
from_spec = mod(ii,Mspecies);
to_vol = ceil(jj/Mspecies);
to_spec = mod(jj,Mspecies);
ok = ok && sum(from_vol ~= to_vol & from_spec ~= to_spec) == 00;
ok = ok && norm(sum(SD,1),inf) < 1e-12;

%-------------------------------------------------------------------------------
function ok = l_spin2
%L_SPIN2 Test of RPARSE_INLINE/SEQEXPAND.

ok = 1;

Nstates = 7;
K0 = reshape(1:Nstates^2,Nstates,Nstates);
umod = rparse_inline([],{'A$i+B$j > K0_$i_$j > C'}, ...
                     {'A$i' 'B$j' 'C'}, ...
                     {'K0_$i_$j' K0}, ...
                     {'i' 1:Nstates 'j' 1:Nstates});

% remove all void reactions
K = umod.inline_propensities.K;
I = umod.inline_propensities.I;
N = umod.N;
ikeep = find(K(1,:));
K = K(:,ikeep);
I = I(:,ikeep);
N = N(:,ikeep);

Nreactions = numel(K0);
[i,j] = ind2sub([Nstates Nstates],1:Nreactions);
ii = [i; j+Nstates; (2*Nstates+1)*ones(1,Nreactions)];
jj = repmat(1:Nreactions,3,1);
ss = repmat([-1 -1 1]',1,Nreactions);
N_ = sparse(ii,jj,ss);

% inline propensities
K_ = [K0(:)'; zeros(2,Nreactions)];
I_ = [i; j+Nstates; ones(1,Nreactions)];
ikeep = find(K_(1,:));
K_ = K_(:,ikeep);
I_ = I_(:,ikeep);
N_ = N_(:,ikeep);

% assertion
ok = ok && all([norm(N-N_,inf) norm(K-K_,inf) norm(I-I_,inf)] < 1e-12);

%-------------------------------------------------------------------------------
