function [param,Theta,dTheta,len] = NDRparam(ncase)
%NDRparams Returns parameters for Notch-Delta-Reporter model.
%   [PARAM,Theta,dTheta,len] = NDRparams(Ncase) returns parameters
%   PARAM for Ncase = 0...8.

% Attempt att reproducing the array of examples in Figure 3/Table 1,
% of "A new mechanism for spatial pattern formation via lateral and
% protrusion-mediated lateral signalling" by Z. Hadjivasiliou,
% G. L. Hunter, and B. Baum in J. R. Soc. Interface
% 13(124):1--10. 2016.

% S. Engblom 2018-01-28 (Revision)
% S. Engblom 2016-12-06

param.betaN = 100;
param.betaD = 500;
param.betaR = 300000;

param.kt = 2;
param.kc = 0.5;
param.kRS = 1e7;

switch ncase
 case 0
  % without protrusions at all
  dTheta = pi;
  len = 0;
  param.wa = 1;
  param.qa = 1;
  param.wb = 0;
  param.qb = 0;
 case 1
  % inspired by figure 3a
  dTheta = pi;
  len = 3.5;
  param.wa = 1;
  param.qa = 1;
  param.wb = 1;
  param.qb = 1;
 case 2
  % inspired by figure 3b
  dTheta = pi;
  len = 3.5;
  param.wa = 1;
  param.qa = 0.001;
  param.wb = 0.06;
  param.qb = 0.06;
 case 3
  % inspired by figure 3c
  dTheta = pi;
  len = 3.5;
  param.wa = 1;
  param.qa = 0.00001;
  param.wb = 0.015;
  param.qb = 0.015;
 case 5
  % inspired by figure 3e
  Theta = pi/2+pi/6;
  dTheta = pi/20;
  len = 5;
  param.wa = 1;
  param.qa = 0.001;
  param.wb = 0.25;
  param.qb = 0.15;
 case 6
  % inspired by figure 3f
  Theta = pi;
  dTheta = pi/20;
  len = 5;
  param.wa = 1;
  param.qa = 0.001;
  param.wb = 0.2;
  param.qb = 0.15;
 case 7
  % inspired by figure 3g
  Theta = -1; % signals perpendicular to radii
  dTheta = pi/20;
  len = 5;
  param.wa = 1;
  param.qa = 0.001;
  param.wb = 0.1;
  param.qb = 0.05;
 case 8
  Theta = -2; % signals radial direction
  dTheta = pi/20;
  len = 5;
  param.wa = 1;
  param.qa = 0.001;
  param.wb = 0.25;
  param.qb = 0.15;
 otherwise
  error('Unknown case.');
end
