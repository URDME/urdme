doffigure = 0;
normalfigure = 0;
deadfigure = 0;

% simulation interval
Tend = 10;
tspan = linspace(0,Tend,101);
% report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

%     The user specified cutoff and rate parameters for the proliferation,
%     death, degradation and consumption rules.

exp = 2;
%Experiments
if exp == 0
    %Normal run------------------------------------
    init = 1;
    start_value = 1;
    radius = 0.1;
    rates_type = 1;
elseif exp ==1
    %initial: relaxation---------------------------
    init = 1;
    start_value = 2;
    radius = 0.2;
    rates_type = 1;
elseif exp==2 
    %initial: simulation---------------------------
    init = 2;
    start_value = 1;
    radius = 0.09;
    rates_type = 2;
end

cutoff_bdof = 0.1;
cutoff_deg = 0.0001;
cutoff_remain = 0.01;

if rates_type == 0
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
elseif rates_type == 1  %relaxation
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0;        % rate of death
    r_degrade = 0;     % rate of degradation for already dead cells
elseif rates_type == 2
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells 
end

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary

% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');

if init == 1
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
elseif init == 2
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
else %all dead
    % initial population: circular blob of dead cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,0,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,0,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,1,[Nvoxels^2 1]);
    U_deadnew = fsparse(ii(:),1,1,[Nvoxels^2 1]); %initiera
end


