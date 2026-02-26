function umod = hes1umod_dlcm(umod)
%Model file for the Hes1 model.
%   UMOD = HES1UMOD(UMOD) constructs (augments) the URDME-struct UMOD
%   with reactions for the Hes1 model.
%
%   See also HES1UMOD2D_RUN.

% E. Blom 2025-08-14 (added $i for dlcm cells)
% S. Engblom 2024-04-10

if nargin == 0, umod = []; end

% rates
rates = hes1_params;

% scale parameters, take into account that hes1_parameters uses units [uM]
avogadro = 6.022e23;
rates.alphan = 1e-6*avogadro*rates.alphan;
rates.KM = 1e-6*avogadro*rates.KM;
rates.Kn = 1e-6*avogadro*rates.Kn;

% species
species = {'D$i' 'N$i' 'M$i' 'P$i' 'n$i'}; % w/out explicit outgoing Dll
% meaning: {Dll, Notch, Hes1 mRNA, Hes1 protein, Ngn2}

% transitions and rates
r = cell(1,10);
r{1} = 'n$i > alphaD*n$i > n$i+D$i';
r{2} = 'D$i > muD*D$i > @';
r{3} = '@ > (sd >= $i)*alphaN/SCALE*ldata[0] > N$i'; % neighbouring signal
r{4} = 'N$i > muN*N$i > @';
r{5} = 'N$i > alphaM*N$i/(1+pow(P$i/(KM*vol),k)) > N$i+M$i';
r{6} = 'M$i > muM*M$i > @';
r{7} = 'M$i > alphaP*M$i > M$i+P$i';
r{8} = 'P$i > muP*P$i > @';
r{9} = '@ > (sd >= $i)*alphan*vol/(1+pow(P$i/(Kn*vol),h)) > n$i';
r{10} = 'n$i > mun*n$i > @';

% note the scaling with SCALE in r{3} above to facilitate bifurcation
% analysis
rates.SCALE = 'gdata'; % passed as global data

% sort it out
umod = rparse(umod,r,species,rates,'hes1.c',{'i', 1:2});

% add some more fields
umod.gdata = 1; % SCALE = 1 (unless changed)
if ~isfield(umod,'vol')
  umod.vol = 1; % scalar unit default
  umod.sd = 1; % a single subdomain
end
umod.u0 = ones(numel(species),numel(umod.vol));
if ~isfield(umod,'tspan')
  umod.tspan = linspace(0,48*60,101);
end
