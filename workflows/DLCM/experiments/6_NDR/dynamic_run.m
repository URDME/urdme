% Notch-delta-reporter model under tissue growth.

% S. Engblom 2018-01-25

% uncomment in blocks to repeat
clear all

% same DLCM-topology used repeatedly below
[P,E,T,gradquotient,V,R,L,dM,N,Nj,Np,U,param] = DLCMlayer(40);

% (1) simulate the growth process first
tstart = tic;
Tend = 200;
growth
time_DLCM = toc(tstart)
save data/dynamic_growth

return;

% (2) ODE
load data/dynamic_growth
tstart = tic;
sol = 1;
NDR_ODEvsSSA_dynamic
time_ODE = toc(tstart)
save data/dynamicODE

return;

% volume variables
VolND = 400;
VolR = 400;

% (3) SSA
load data/dynamic_growth
sol = 2;

% case 1 (default)
tstart = tic;
NDR_ODEvsSSA_dynamic
time_SSA = toc(tstart)
save data/dynamicSSA1

% case 2
param = NDRparam(2);
NDR_ODEvsSSA_dynamic
save data/dynamicSSA2

% case 6
[param,Theta,dTheta,len] = NDRparam(6);
% polarized protrusion in this case:
[P,E,T,gradquotient,V,R,L,dM,N,Nj,Np,U,~] = DLCMlayer(40,Theta,dTheta,len);
NDR_ODEvsSSA_dynamic
save data/dynamicSSA6

return;

% (4) RDME

% case 1 (default)
load data/dynamic_growth
tstart = tic;
NDR_RDME_dynamic
time_RDME = toc(tstart)
save data/dynamicRDME1

% case 2
param = NDRparam(2);
tstart = tic;
NDR_RDME_dynamic
time_RDME = toc(tstart)
save data/dynamicRDME2

% case 6
[param,Theta,dTheta,len] = NDRparam(6);
% polarized protrusion in this case:
[P,E,T,gradquotient,V,R,L,dM,N,Nj,Np,U,~] = DLCMlayer(40,Theta,dTheta,len);
NDR_RDME_dynamic
save data/dynamicRDME6

return;
