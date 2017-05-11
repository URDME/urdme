function umod = urdme2comsol(umod,verbose)
%URDME2COMSOL Update .comsol-field with URDME result.
%   UMOD = URDME2COMSOL(UMOD,VERBOSE) inserts the trajectory
%   (UMOD.U,UMOD.tspan) into UMOD.comsol under the verbosity VERBOSE.
%
%   The solution matrix U produced by URDME contains integers (number
%   of species), not concentrations. URDME2COMSOL turns these values
%   into concentration for handling by the Comsol API. On return
%   UMOD.comsol.sol1 is the solution thus added which corresponds to
%   input (UMOD.U,UMOD.tspan).
%
%   See also URDME, COMSOL2URDME.

% S. Engblom 2017-02-15 (Major revision, URDME 1.3, Comsol 5)
% P.Bauer 2012-08-30 (urdme2comsol)
% V.Gerdin 2012-02-27 (rdme2mod)
% A. Hellander 2010-05-04
% J. Cullhed 2008-08-06 (rdme2fem)

if nargin < 2, verbose = 0; end

% extend mesh from model
l_info(verbose,1,'Start mphxmeshinfo...');
xmi = mphxmeshinfo(umod.comsol);
l_info(verbose,1,' ...done.\n');

% number of species and number of cells
nodes = xmi.nodes;
[Mspecies,Ncells] = size(nodes.dofs);
Ndofs = Mspecies*Ncells;

% permutation of species
if isfield(umod,'species')
  fieldnames = umod.species;
else
  fieldnames = xmi.fieldnames;
end
modName = char(umod.comsol.modelNode.tags());
if isempty(strfind(char(fieldnames(1)),modName))
  fieldnames = strcat(modName,'.',fieldnames);
end
[foo,sp] = ismember(fieldnames,nodes.dofnames);
isp(sp) = 1:numel(sp); % invert order

% change from absolute figures to concentration
U = reshape(umod.U,Mspecies,Ncells,[]);
V = repmat(umod.vol(:)'*6.022e23,[Mspecies 1]);
for i = 1:numel(umod.tspan)
  U(:,:,i) = U(:,:,i)./V;
end
U = reshape(U(isp,:,:),Ndofs,[]);

% store the resulting U in umod.comsol
l_info(verbose,1,'Saving U to model...');
for i = 1:numel(umod.tspan)
  umod.comsol.sol('sol1').setU(i,U(:,i));
end
umod.comsol.sol('sol1').setPVals(umod.tspan);
umod.comsol.sol('sol1').createSolution;
l_info(verbose,1,' ...done.\n');

%-------------------------------------------------------------------------
function l_info(level,lim,msg)
%L_INFO Display information.
%   L_INFO(LEVEL,LIM,MSG) Displays the message MSG whenever LEVEL >=
%   LIM.

if level >= lim, fprintf(msg); end

%-------------------------------------------------------------------------
