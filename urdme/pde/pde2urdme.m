function umod = pde2urdme(P,T,Dexpr,verbose)
%PDE2URDME Create URDME-structure from a PDE Toolbox model.
%   UMOD = PDE2URDME(P,T,Dexpr,VERBOSE) creates an URDME-structure
%   UMOD from the triangle mesh (P,T) and diffusion expression Dexpr
%   under the verbosity VERBOSE.
%
%   Dexpr = {'D_spec1' 'Dspec2' ...} is a cell-vector containing
%   diffusion expressions for the different species. The diffusion
%   expression may involve x, y, and sd (subdomain number).
%
%   Cautionary: use the same ordering in Dexpr as in the specification
%   of the stoichiometric matrix UMOD.N.
%
%   Note that UMOD.sd contains an interpolated value of the subdomain
%   number. You will need to modify this value to arrive at an integer
%   value.
%
%   See also URDME, URDME2PDE.

% S. Engblom 2017-02-21

if nargin < 4, verbose = 0; end

% number of species and number of cells
Mspecies = numel(Dexpr);
Ncells = size(P,2);
Ndofs = Ncells*Mspecies;

% assemble matrices
l_info(verbose,1,'Starting assembly...');
D = cell(size(Dexpr));
% evaluate element volume and subdomain first
[~,M,SD] = assema(P,T,0,1,'sd');
vol = full(sum(M,1))';

% assemble all species
diagsum1 = 0;
I = zeros(0,1); J = zeros(0,1); S = zeros(0,1);
for i = 1:Mspecies
  D{i} = assema(P,T,Dexpr{i},0,0);

  % explicitly invert the lumped mass matrix and filter the diffusion matrix
  [ii,jj,ss] = find(D{i});
  diagsum1 = diagsum1+sum(diag(D{i})./vol);
  ss = -ss./vol(jj);
  ixkeep = find(ss > 0);
  I = [I; (ii(ixkeep)-1)*Mspecies+i];
  J = [J; (jj(ixkeep)-1)*Mspecies+i];
  S = [S; ss(ixkeep)];
end
l_info(verbose,1,' ...done.\n');

% rebuild diffusion matrix
D = sparse(I(:),J(:),S(:),Ndofs,Ndofs);
d = full(sum(D,1));
D = D+sparse(1:Ndofs,1:Ndofs,-d);
diagsum2 = sum(d);

% check if the difference is too big
if abs(diagsum2-diagsum1) > abs(diagsum1)*0.10
  warning('Many off diagonal negative elements in D.');
end

% this creates interpolated values of the subdomain number
sd = SD./vol;
% (typically to be fixed to arrive at integer values)

% finalize the PDE-part of the URDME structure
umod = struct('D',D, ...
              'vol',vol, ...
              'sd',sd, ...
              'pde',struct('P',P,'T',T,'Dexpr',{Dexpr}));

%-------------------------------------------------------------------------
function l_info(level,lim,msg)
%L_INFO Display information.
%   L_INFO(LEVEL,LIM,MSG) Displays the message MSG whenever LEVEL >=
%   LIM.

if level >= lim, fprintf(msg); end

%-------------------------------------------------------------------------
