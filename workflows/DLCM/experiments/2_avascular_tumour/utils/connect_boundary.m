function [extdof] = connect_boundary(extdof, M, N, P)
% CONNECT_BOUNDARY Connects FEM boundary nodes
% Connects the boundary nodes of a subset of nodes in a regular, 
% triangular mesh by including the missing nodes. The boundary nodes
% are assumed to comprise a circle of radius R.
% This function adds the 'furthest' nodes such that area of the domain
% enclosed by the boundary nodes is maximized, but minimizing the effect 
% of the skewed boundary for Poisson equations solved within the domain.
% 
%   input:  extdof  -   Boundary nodes, found externally as the nodes closest 
%                       to R
%           M       -   Mass matrix
%           N       -   Neighbor matrix (Adjacency matrix + identity matrix)
%           P       -   Vector of node coordinates, by Matlabs P, E, T standard
% 
%   output: extdof     - new set of connected boundary nodes
%

% E. Blom 2022-09-01

R = 1; % radius of boundary nodes.
extvec = full(fsparse(extdof,1,1 ,[size(M,1) 1]));  % extdof boolean vector

% fill in boundary connectivity gaps
% find boundary edges
bdof = (N*extvec <= 2 & extvec); % alone or endpoint extdofs...
cdof = find(N*bdof == 2 & ~extvec); % ... nodes connecting these
tbdof = ((N*bdof).*bdof == 2);  % troublesome connected bdofs...
tcdof = find(N*tbdof == 2 & ~extvec); % ...new nodes that connect these

% find the nodes that SHOULD connect 2 bdofs, too keep
% (Quite complicated approach - might be able to simplify!)
tN = N(tbdof,tcdof);        % local neighbour-matrix for troublesome dofs
tii = find(sum(tN,2)==4);   % index for 'false' troublesome dofs...
[~, b] = find(tN(tii,:));   % ... and their troublesome neighbours 
[uniqueb, i, j] = unique(b,'first');    % find overlapping neighbors...
idx_dupl = find(not(ismember(1:numel(b),i))); %...i.e., duplicates in b
cdof_keep = tcdof(b(idx_dupl)); % keep these nodes...
tcdof = setdiff(tcdof, cdof_keep); %...i.e., remove these from list of removals

cdof = setdiff(cdof, tcdof);    % remove these indices

%indof = cdof((P(1,cdof).^2 + P(2,cdof).^2) < 1);
cdof = cdof((P(1,cdof).^2 + P(2,cdof).^2) > R);
% adding all cdof to boundary gives smallest connected boundary,
% but note that there likely exists 'bundled' boundary dofs with 
% > 2 connections to other boundary dofs!
extdof = union(extdof, cdof);
end

