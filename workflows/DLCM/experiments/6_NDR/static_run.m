% Notch-delta-reporter model: ODE, SSA, RDME.
%
%   This file runs 3 models for calibration purposes.

% S. Engblom 2018-01-24

% uncomment in blocks to repeat
clear all

% $$$ % same DLCM-topology used throughout
% $$$ [P,E,T,gradquotient,V,R,L,dM,N,Nj,Np,U,param] = DLCMlayer(20);
% $$$ save data/staticDLCM
% $$$ 
% $$$ return;

% volume variables
VolND = 400;
VolR = 400;

% (1) ODE
% $$$ load data/staticDLCM
% $$$ sol = 1;
% $$$ Tend = 200;
% $$$ NDR_ODEvsSSA_static
% $$$ save data/staticODE
% $$$ 
% $$$ return;

% (2) SSA
% $$$ load data/staticDLCM
% $$$ sol_restart = false;
% $$$ sol = 2;
% $$$ Tend = 200;
% $$$ NDR_ODEvsSSA_static
% $$$ save data/staticSSA
% $$$ 
% $$$ return;

% (3) RDME
% $$$ load data/staticDLCM
% $$$ % due to the initial transient, this is extremely slow:
% $$$ sol_restart = false;
% $$$ % but this is faster as it avoids this transient:
% $$$ sol_restart = true;
% $$$ Tend = 200;
% $$$ NDR_RDME_static
% $$$ save data/staticRDME
% $$$ 
% $$$ return;

% plot dynamics of selected voxel
figure(1), h = gca;
while 1
  % fetch voxel
  axes(h);
  [x,y,button] = ginput(1);
  if button == 3, break; end
  [~,ii] = min((P(1,:)-x).^2+(P(2,:)-y).^2);
  ii

  % plot the dynamics of that voxel
  NN = zeros(1,numel(tspan));
  DD = zeros(1,numel(tspan));
  RR = zeros(1,numel(tspan));
  for kk = 1:numel(tspan)
    NN(kk) = Nsave{kk}(ii,1);
    DD(kk) = Dsave{kk}(ii,1);
    RR(kk) = Rsave{kk}(ii,1);
  end 
  figure(4), clf,
  subplot(3,1,1);
  plot(tspan,NN,'bo-'); title('N'); ax1 = ylim;
  subplot(3,1,2);
  plot(tspan,DD,'ro-'); title('D'); ax2 = ylim;
  subplot(3,1,3);
  plot(tspan,RR,'go-'); title('R'); ax3 = ylim;
end
