function umod = hes1umod(umod)
%Model file for the Hes1 model.
%   UMOD = HES1UMOD(UMOD) constructs (augments) the URDME-struct UMOD
%   with reactions for the Hes1 model.
%
%   See also HES1UMOD2D_RUN.

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
species = {'D' 'Din' 'N' 'M' 'P' 'n'};
% meaning: {Dll, outgoing Dll, Notch, Hes1 mRNA, Hes1 protein, Ngn2}

% transitions and rates
r = cell(1,10);
r{1} = 'n > alphaD*n > n+D';
r{2} = 'D > muD*D > @';
% creation of N from D at rate alphaN is split into (1) a creation of
% Din followed by (2) a diffusion event of a Din into an N in the
% neighbor voxel
r{3} = 'D > alphaN/SCALE*D > D+Din';
r{4} = 'N > muN*N > @';
r{5} = 'N > alphaM*N/(1+pow(P/(KM*vol),k)) > N+M';
r{6} = 'M > muM*M > @';
r{7} = 'M > alphaP*M > M+P';
r{8} = 'P > muP*P > @';
r{9} = '@ > alphan*vol/(1+pow(P/(Kn*vol),h)) > n';
r{10} = 'n > mun*n > @';

% note the scaling with SCALE in r{3} above to facilitate bifurcation
% analysis
rates.SCALE = 'gdata'; % passed as global data

% sort it out
umod = rparse(umod,r,species,rates,'hes1.c');

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
