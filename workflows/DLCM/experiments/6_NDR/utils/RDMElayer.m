function [P,E,T,G] = RDMElayer
%RDMElayer Utility function which sets up the RDME-layer.
%   [P,E,T,G] = RDMElayer creates a crude single cell geometry G
%   discretized with triangulation (P,E,T) for use by PDE2URDME.

% S. Engblom 2018-01-24

% create the geometry, composed of two circles
C1 = [1 0 0 1]';
C2 = [1 0 0 0.45]';
G = decsg([C1 C2],'C1+C2',char('C1','C2')');

% create the mesh
[P,E,T] = initmesh(G,'hmax',0.5);
