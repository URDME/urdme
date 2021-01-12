% Simulation of an avascular tumour model.
%
%   Avascular tumour growth: An initial circular population cells (one
%   per voxel) lie in a domain rich in oxygen. Cells consume oxygen at
%   a constant rate, lambda. Cells occupying a voxel with oxygen above
%   cutoff_prol can proliferate at a rate r_prol. Cells occupying
%   voxels with an oxygen concentration below cutoff_die can die at a
%   rate r_die.  Dead cells are represented with a voxel with value
%   -1, these dead cells can degrade and stop occupying space at a
%   rate r_degrade.
%
%   Permeability: Drate1 describes the rate diffusion rate of tumour
%   cells invading previously unvisited voxels. Drate2 is the rate
%   cells move into previously occupied but currently empty
%   voxels. Drate3 is the rate cells move into voxels that are already
%   occupied.

% C. Jayaweera & A. Graf Brolund 2020-12(revision)
% S. Engblom 2017-12-27 (revision)
% D. B. Wilson 2017-09-05
% S. Engblom 2017-02-11

clear;
clc;
close all;

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered
%Nvoxels = 243; 

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

D=1; %D_rate 

% simulation interval
Tend = 100;
tspan = linspace(0,Tend,101);
timescaling = 0.005;

% initial population: circular blob of living cells    
start_value = 10;
radius = 0.07;
r = sqrt(P(1,:).^2+P(2,:).^2);
ii = find(r < radius); % radius of the initial blob
U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize

% parameters
cons = 0.0015;        % consumption of oxygen by cells
cutoff_bdof = 0.1;      % lower bound for bdof

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T);       %N gives the neighbours 
neigh = full(sum(N,2));

% dofs for the sources at the extreme outer circular boundary
[xc,yc] = getmidpointcircle(1/2*(Nvoxels+1),1/2*(Nvoxels+1),1/2*(Nvoxels-1));
irem = find(xc < 1 | yc < 1 | xc > Nvoxels | yc > Nvoxels);
xc(irem) = [];
yc(irem) = [];
extdof = find(sparse(xc,yc,1,Nvoxels,Nvoxels));

% visit marker matrix: 1 for voxels who have been occupied
VU = (U ~= 0);

% representation of solution: cell-vector of sparse matrices
Usave = cell(1,numel(tspan));
Usave{1} = U;
Udsave = cell(1,numel(tspan));
Udsave{1} = U_dead;

Oxysave = cell(1,numel(tspan));
bdofsave = cell(1,numel(tspan));
sdofsave = cell(1,numel(tspan));
sdofbsave = cell(1,numel(tspan));

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

while tt <= tspan(end)
    U = U_new;
    U_dead = U_deadnew;
    
    dof_calculation;           %Calculation of dofs
    PrOx_calculation;          %Laplacian calculation

    move_calculations;         %sdof and bdof movement calculations
 
    ind_rates_sdof_n = find(rates_sdof(sdof_m_)<0);
    ind_rates_bdof_n = find(rates_bdof(bdof_m_)<0);

    dt_sdof = U_new(sdof_m(ind_rates_sdof_n))./(-rates_sdof(sdof_m_(ind_rates_sdof_n)));
    dt_bdof = U_new(bdof_m(ind_rates_bdof_n))./(-rates_bdof(bdof_m_(ind_rates_bdof_n)));

    dt = min([dt_sdof; dt_bdof;(0.1*Tend)])*timescaling;
   
    tspan_calculation;         %save times series of current states
    
    %Euler step
    
    % sdof_m
    U_new(Adof) = U_new(Adof) + rates_sdof.*dt;

    % bdof_m
    U_new(Adof) = U_new(Adof) + rates_bdof*dt;
    
    tt = tt+dt;
    report(tt,U,'');
    
    % update the visited sites
    VU = VU | U;
end
report(tt,U,'done');

%%
TumorGraphics;
