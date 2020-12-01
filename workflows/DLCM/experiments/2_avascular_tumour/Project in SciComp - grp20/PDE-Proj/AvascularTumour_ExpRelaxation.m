% Avascular tumour model
% Experiment: Avascular 
% Chrishani Jayaweera 2020-12-01 

clear;
clc;
close all;

%profile on

init_experiment;

% assemble minus the Laplacian on this grid (ignoring BCs), the voxel
% volume vector, and the sparse neighbor matrix
[L,dM,N] = dt_operators(P,T); %N=raden ger alla voxlar, 1or fÃ¶r de som Ã¤r grannar med radens voxel, tom eller ej
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

birth_count = 0;
tt = tspan(1);
i = 1;
% logic for reuse of LU-factorizations
updLU = true;
La = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
OLa = struct('X',0,'L',0,'U',0,'p',0,'q',0,'R',0);
% event counter
% Ne = struct('moveb',0,'moves',0,'birth',0,'death',0,'degrade',0);

% oxygen Laplacian
OLa.X = L;
OLai = fsparse(extdof,extdof,1,size(OLa.X));
OLa.X = OLa.X-OLai*OLa.X+OLai;   %add Dirichlet(?) oxygen at the oxygen source/outer circle
[OLa.L,OLa.U,OLa.p,OLa.q,OLa.R] = lu(OLa.X,'vector');

while tt <= tspan(end)
    U = U_new;
    U_dead = U_deadnew;
    
    dof_calculation;           %Calculation of dofs

    RHS_calculation;           %Laplacian calculation

    intensity_calculation;    %Calculate intensties of rates of events
    lambda = sum(intens);
    dt = 1/lambda;
    
    tspan_calculation; %save snapshots to time vector
    
    move_calculations; %sdof and bdof calculations
    
    change_calculation;  %Proliferation, death and degradation
    
    updLU = true; % boundary has changed
   
    tt = tt+dt;
    report(tt,U,'');
    
    % update the visited sites
    %VU = VU | U;
end
report(tt,U,'done');

% return;
% profile off
% profile report

%%
TumorGraphics;

