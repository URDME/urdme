function dlcm
%DLCM Discrete Laplacian Cell Mechanics solver.
%   The URDME DLCM solver is an implementation of the cell population
%   framework described in [1], introduced in [2]. Cells reside in
%   voxels on a spatial grid, their states updated according to
%   population-level reaction-transport processes based on their local
%   micro-environment. Indepedent reaction processes may also occur in
%   each cell using URDME's SSA solver internally.
%
%   The population-level states and events are defined directly in
%   the umod-struct, as are the micro-environment quantity models
%   (see example below). These quantities must be defined as Q# for
%   quantitiy number #, starting from Q1; all reaction events at this
%   level must be defined before the quantity models. This ensures
%   correct parsing through DLCM2URDME.
%
%   In addition to the required URDME fields, the solver requires
%   solver arguments as detailed below, see DLCM2URDME for how to
%   conveniently construct these. The formally required URDME field G
%   is not used by the solver while D is passed as solver argument.
%
%   Required fields
%   -----------------------------------------------------------------------
%   P            Voxel position (Ndim-by-Nvoxels matrix)
%   gradquotient Ratio between shared edge length to distance between
%                neighboring voxels (scalar or Nvoxels-by-Nvoxels sparse
%                matrix)
%   Here Ndim is the number of spatial dimensions and Nvoxels is the number
%   of voxels.
%
%   Required fields with defaults
%   -----------------------------------------------------------------------
%   D            Discrete Laplacian
%   Ne           Mesh neighbor matrix
%   Rates        Migration rates (function returning 1-by-Nmig cell vector)
%   Drate        Scaling of the migration rates array function in Rates 
%                (function returning 1-by-Nmig cell vector)
%   Here Nmig is the number of migration potentials used in the model.
%
%   Optional fields
%   -----------------------------------------------------------------------
%   T                 Mesh triangulation connectivity list
%                     (m-by-n array, see INITMESH or meshToPet)
%   sigma             Surface tension coefficients (symmetric
%                     Ntypes+1-by-Ntypes+1 matrix)
%   mumod             URDME umod struct for internal states and events
%   maxdt_fun         Maximum SSA time step
%   ldata_fun         Internal event rates function (function returning
%                     1-by-Nldata cell vector)
%   internal_state    Flag for continuous or discrete (default) internal
%                     states ({'discr'} | 'cont')
%   private.dlcm.glix Global cell indices, only output 
%                     (2-by-Nvoxels-by-Ntime)
%
%   Here Ntypes is the number of cell types in the model and Nldata is
%   the number of complex internal rate functions (regarded by mumod
%   as ldata). Ntime is the number of time stamps where simulation
%   data is saved.
%
%   Examples:
%     %% (1) geometry
%     Nvoxels = 41;
%     mesh_type = 1;  % cartesian mesh
%
%     % Simulate to Tend and save states at Tres intervals
%     Tend = 10000;
%     Tres = 100;
%     ntypes = 2; % number of cell types: here two arbitrary types
%
%     [P,E,T,gradquotient] = basic_mesh(mesh_type,Nvoxels);
%     [V,R] = mesh2dual(P,E,T,'voronoi'); % dual mesh for visualizaton
%
%     % Get boundary dofs, extdof, where the micro-environment BCs are
%     % defined
%     hmax = 2/(Nvoxels-1);               % half voxel size rough estimate
%     extdof = find(P(1,:).^2+P(2,:).^2 > (1-hmax)^2 & ...
%                   P(1,:).^2+P(2,:).^2 <= 1^2);
%
%     %% (2) Migration rates per cell
%     Rates = @(U,Q,QI,P,t){Q(:,1)};   % here, only migration pressure
%     Drate = @(Uf,Ut,Q,QI,P,t){1.*(Uf==1).*(Ut==0)+1.*(Uf==2).*(Ut<2)};
%     % Note: solver considers only migration from
%     % bdof_m and sdof_m, see MEXDLCM
%
%     %% (3) Form population
%     % initial small population
%     ii1 =  find(abs(P(1,:)) < 0.4 & abs(P(2,:)) <= 0.4);  % cell type 1
%     ii2 = [];                                             % cell type 2
%     % U is Ntype-by-Ncells sparse vector, representing the cell population
%     U(1,:) = fsparse(ii1(:),1,2,[Nvoxels^2 1]); % doubly occupied, type 1
%     U(2,:) = fsparse(ii2(:),1,1,[Nvoxels^2 1]); % initially none, type 2
%
%     %% (4) "outer" URDME-struct
%     % Defining reaction events and quantity models here
%     nquants = 1; % number of field states pressure and nutrient
%     Dexpr = cell(1,nquants+ntypes);
%     Dexpr(:) = {1};
%     umod = pde2urdme(P,T,Dexpr);              % construct D matrix, etc.
%     % 'UL' means cell of phenotype L.
%     % Define all reaction events first, and after that quantities:
%     umod = rparse(umod, ...
%               % Define reaction events first
%               {'U1 > 0.1 > U1+U1' ...         % birth event
%                'U1 > 0.05 > U2' ...           % switch phenotype
%                'U2 > 0.01 > U1+U1' ...        % birth+switch
%                'U2 > 0.1 > @' ...             % death event
%               % Define quantity models last
%                'Q1 > (U1>1) > Q1+Q1'}, ...    % Q1 sources = (U1>1)
%               {'U1' 'U2' 'Q1'}, ...
%               {}, ...
%               'test');
%     umod.u0 = [full(U); zeros(1,Nvoxels^2)];
%     umod.sd = ones(1,Nvoxels^2);
%     umod.sd(extdof) = 0;                      % sd encodes boundary dofs
%     umod.tspan = linspace(0,Tend,Tres);       % time steps
%     % load essentials
%     umod = dlcm2urdme(umod,P,gradquotient,[],[],[], ...
%                       'Rates',Rates,'Drate',Drate);
%
%     %% (6) solve
%     umod = urdme(umod,'solver','dlcm','solve',1,'seed',123);
%
%   See also URDME, SSA, DLCM2URDME.
%
%   References:
%     [1] E. Blom, S. Engblom. "DLCM: a versatile multi-level solver for
%     heterogeneous multicellular systems". ArXiv Preprint (2025)
%     [2] S. Engblom, D. B. Wilson, and R. E. Baker.
%     "Scalable population-level modelling of biological cells
%     incorporating mechanics and kinetics in continuous time."
%     Royal Society open science 5.8 (2018).
%

% E. Blom 2024-04-25

error('This file should not be called directly.');
