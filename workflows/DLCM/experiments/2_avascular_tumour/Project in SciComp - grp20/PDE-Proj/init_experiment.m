
% cells live in a square of Nvoxels-by-Nvoxels
Nvoxels = 121; % odd so the BC for oxygen can by centered

% fetch Cartesian discretization
[P,E,T,gradquotient] = basic_mesh(1,Nvoxels);  %gradquotient=1 for Cartesian mesh
[V,R] = mesh2dual(P,E,T,'voronoi');


D=1; %D_rate 

% simulation interval
Tend = 90;
tspan = linspace(0,Tend,101);
% report(tspan,'timeleft','init'); % (this estimator gets seriously confused!)

%     The user specified cutoff and rate parameters for the proliferation,
%     death, degradation and consumption rules.
exp = 0;

%Experiments
if exp == 0
    % Normal run------------------------------------
    start_value = 1;
    radius = 0.09;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells
    
    % Initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]);      %Intialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]);   %Initialize
    
elseif exp ==1
    % Initial: relaxation---------------------------
    start_value = 10;
    radius = 0.07;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0;        % rate of death
    r_degrade = 0;     % rate of degradation for already dead cells
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
elseif exp==2 
    %initial: simulation---------------------------
    start_value = 1;
    radius = 0.09;
    
    cons = 0.0015;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.125;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.125;        % rate of death
    r_degrade = 0.01;     % rate of degradation for already dead cells 
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initialize
    
elseif exp==3 
    %initial: Quick, concentrated------------------
    start_value = 1;
    radius = 0.09;
    
    %Tend=30;
    cons = 0.15;        % consumption of oxygen by cells
    cutoff_prol = 0.65;   % the minimum amount of oxygen for proliferation
    r_prol = 0.3;       % rate of proliferation of singly occupied voxels
    cutoff_die = 0.55;    % the maximum amount of oxygen where cells can die
    r_die = 0.3;        % rate of death
    r_degrade = 0.05;     % rate of degradation for already dead cells 
    
    % initial population: circular blob of living cells
    r = sqrt(P(1,:).^2+P(2,:).^2);
    ii = find(r < radius); % radius of the initial blob
    U = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_new = fsparse(ii(:),1,start_value,[Nvoxels^2 1]);
    U_dead = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
    U_deadnew = fsparse(ii(:),1,0,[Nvoxels^2 1]); %initiera
end

cutoff_bdof = 0.1;
cutoff_deg = 0.0001;
cutoff_remain = 0.01;

% boundary conditions
OBC1 = 0; % BC for the oxygen equation for unvisited boundary
OBC2 = 0; % BC for the visited boundary


% Other experiments
% %All dead initiliazation
%     % initial population: circular blob of dead cells
%     r = sqrt(P(1,:).^2+P(2,:).^2);
%     ii = find(r < radius); % radius of the initial blob
%     U = fsparse(ii(:),1,0,[Nvoxels^2 1]);
%     U_new = fsparse(ii(:),1,0,[Nvoxels^2 1]);
%     U_dead = fsparse(ii(:),1,1,[Nvoxels^2 1]);
%     U_deadnew = fsparse(ii(:),1,1,[Nvoxels^2 1]); %initialize



