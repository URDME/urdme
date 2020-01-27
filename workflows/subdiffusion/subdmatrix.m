function SD = subdmatrix(A,gamma,D,select)
%SUBDMATRIX Compute generalized diffusion matrix.
%   SD = SUBDMATRIX(A,GAMMA,D) returns a generalized diffusion matrix
%   SD for A an internal states transfer matrix, GAMMA the scalar
%   diffusion constants for the internal states, and D the diffusion
%   operator on the mesh considered.
%
%   SD = SUBDMATRIX(A,GAMMA,D,SELECT) additionally selects the species
%   to expand into internal states through the logical vector SELECT;
%   SELECT(i) is true if species i is to be expanded and false
%   otherwise. The total number of species before expanding is hence
%   NUMEL(SELECT).
%
%   No error-checking is performed.

% S. Engblom 2017-04-17 (added select)
% S. Engblom 2017-03-21 (coarse-grained version)
% S. Engblom 2014-05-09 (getdmatrix)

if nargin < 4
  Nstates = numel(gamma);
  Nmesh = size(D,1);

  % eq. (2.26) \ref{eq:diffequ1}: combine the two to produce a single
  % generalized diffusion matrix
  SD = kron(D,sparse(1:Nstates,1:Nstates,gamma(:)'))+ ...
       kron(sparse(1:Nmesh,1:Nmesh,1),A);
  % (SD will be such that sum(SD,1) is 0 provided (A,D) does this)
else
  select = logical(select);
  Mspecies = numel(select);
  Nstates = numel(gamma);
  Nmesh = size(D,1)/Mspecies;

  % separate into untouched part and internal states part
  ix = reshape(1:Mspecies*Nmesh,Mspecies,Nmesh);
  jx = ix(~select,:);
  ix = ix(select,:);
  % untouched part
  D1 = D(jx,jx);
  % internal states part
  D2 = kron(D(ix,ix),sparse(1:Nstates,1:Nstates,gamma(:)'))+ ...
       kron(sparse(1:Nmesh*sum(select),1:Nmesh*sum(select),1),A);
  
  % construct the associated permutation of all species, including
  % internal states
  ix = 1:Nstates:Mspecies*Nstates;
  jx = ix(~select);
  ix = ix(select);
  ix = repmat(ix,Nstates,1);
  [~,p] = sort([jx(:); ix(:)]); % note: sort is stable

  % the total permutation then follows, albeit nontrivially:
  ix = [reshape(1:size(D1,1),[],Nmesh); ...
        reshape((1:size(D2,1))+size(D1,1),[],Nmesh)];
  p = ix(p,:);

  % concatenate and permute accordingly
  SD = blkdiag(D1,D2);
  SD = SD(p,p);
end
