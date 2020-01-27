function channels = squid_membrane(model)
%SQUID_MEMBRANE
%   CHANNELS = SQUID_MEMBRANE(MODEL) creates the struct CHANNELS
%   required for cable_description of the neuron. The input model is
%   to contain all the necessary parameters for the HH neuron, see
%   GIAN_SQUID_STRUCT below.

% A. Senek 2017-05-31

params = giant_squid_struct();

nVoxels = model.nVoxels;
dA = model.dA;  % Adjacency matrix
L = model.L;    % Segment lengths
D = model.D;    % Segment diameters
surf_area = model.surfacearea;  % Segment surface-areas

%% Cytoplasm
channels.R_cytoplasm = params.Neuron.Cytoplasm.Resistivity;
channels.G_axial = gaxial(L,D(2:end),dA(2:end,2:end),channels.R_cytoplasm);

%% Membrane model
channels.C_m = params.Neuron.Membrane.SpecificCapacitance.*surf_area;
channels.Erest = params.Neuron.Membrane.RestingPotential.*ones(nVoxels,1);

channels.G_m = params.Neuron.Membrane.SpecificConductivity.*surf_area;
%% Voltage-gated ion channel models
channels.Erev(1,1) = params.Ionchannel.K.ReversalPotential;
channels.density(1,1) = params.Ionchannel.K.ChannelDensity;
channels.G_single(1,1) = params.Ionchannel.K.SingleChannelConductance;

channels.Erev(2,1) = params.Ionchannel.Na.ReversalPotential;
channels.density(1,2) = params.Ionchannel.Na.ChannelDensity;
channels.G_single(2,1) = params.Ionchannel.Na.SingleChannelConductance;

%% Leakage ion channel model
channels.lchannel.G = params.LeakChannel.SpecificConductivity.*surf_area;
channels.lchannel.Erev = params.LeakChannel.ReversalPotential.*ones(nVoxels,1);

%-------------------------------------------------------------------------
function Ga = gaxial(L, D, dA, resistivity)
%   Ga = GAXIAL(L,D,dA,resistivity)
% INPUT:
%   L  (1,nVoxels)             length of voxels
%   D  (1,nVoxels)            diameters of voxels (D of tree minus root)
%   dA (nVoxels,nVoxels)   directed adjacency matrix of voxels
%
%   OUTPUT:
%   Ga(nVoxels,nVoxels), contains the conductance between neuronal compartments

% Index to direct parent:
parent = dA*(1:size(dA,2))';

% Assume that the currents are parallel to the axis and homogenous,
% approximate the conductances as in series going from the middle of the voxels
% g = 1/(rho * [1/2 L/A + 1/2 l/a])
N = numel(L);
counter = 0;
i = []; j = []; ga = [];
for s=1:N %voxel numbers
  Ls = L(s);
  Ds = D(s);
  parents = parent(s);
  if(parents > 0)
    counter = counter+1;
    i(counter) = s;
    j(counter) = parents;
    ga(counter) = 1/(resistivity/2* (Ls/(pi*(Ds/2)^2) + L(parents)/(pi*(D(parents)/2)^2)) );
  end
  for child=find(dA(:,s)==1)'
    counter = counter+1;
    i(counter) = s;
    j(counter) = child;
    ga(counter) = 1/(resistivity/2 * (Ls/(pi*(Ds/2)^2) + L(child)/(pi*(D(child)/2)^2)) );
  end
end
Ga = sparse(i,j,ga,N,N);

%-------------------------------------------------------------------------
function giant_squid = giant_squid_struct()
% giant_squid = GIANT_SQUID_STRUCT()
%
%   Returns a struct which contains all the relevant parameters regarding
%   the giant squid neuron. For reference see Rallpack 3 at:
%   https://github.com/borismarin/genesis2.4gamma/tree/master/rallpack

giant_squid.Neuron.Cytoplasm.Resistivity = 100 * 1e-2; %   Resistivity (ohm) (centi meter)

giant_squid.Neuron.Membrane.SpecificCapacitance = 1e-5;%   (micro farad) (centi meter)^(-2)
giant_squid.Neuron.Membrane.RestingPotential = -65;% (milli volt)
giant_squid.Neuron.Membrane.SpecificConductivity = 0;% (siemens) (centi meter)^(-2)

giant_squid.Ionchannel.Na.ChannelDensity = 330;% (micro meter)^(-2)
giant_squid.Ionchannel.Na.ReversalPotential = 50;% (micro meter)^(-2)
giant_squid.Ionchannel.Na.SingleChannelConductance = 1200/330 * 1e-6;% (pico siemens)

giant_squid.Ionchannel.K.ChannelDensity = 30;% (micro meter)^(-2)
giant_squid.Ionchannel.K.ReversalPotential = -78;% (micro meter)^(-2)
giant_squid.Ionchannel.K.SingleChannelConductance = 12 * 1e-6;% (pico siemens)

giant_squid.LeakChannel.SpecificConductivity = 1/40000 * 1e-2;% (siemens) (centi meter)^(-2)
giant_squid.LeakChannel.ReversalPotential = -65;% (milli volt)

%-------------------------------------------------------------------------
