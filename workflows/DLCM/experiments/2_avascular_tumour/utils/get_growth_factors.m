function [Lambda, sigma_stable] = get_growth_factors(args)
% GET_GROWTH_FACTORS evaluates the analytical 
%   growth factors of each perturbation mode in modes
%
%   in: args.rp - radius of tumor per time step
%       args.rq - radius of quiescent-proliferating interface per time step
%       args.rn - radius of per time step
%       args.sigma - surface tension strength
%       args.mu_prol - rate of proliferation
%       args.mu_die - rate of death
%       args.modes - integer perturbation modes to be evaluated
%
%   out: Lambda - amplification factors for each input perturbation mode
%        sigma_stable - the sigma required to stabilise mode k (stationary state)

% E. Blom 2023-03-09

% 'load' parameters
rp = args.rp; rq = args.rq; rn = args.rn;
sigma = args.sigma; mu_prol = args.mu_prol; mu_die = args.mu_death;
modes = args.modes;

% get lambda vs rp or lambda vs sigma
if numel(rp) > 1
Lambda = zeros(numel(modes), numel(rp));
elseif numel(sigma) > 1
Lambda = zeros(numel(modes), numel(sigma));
end

for k = modes
    % get estimated regional perturbations
    rqc = -1./k*rp.*(rp.^-k - rp.^k)./(rn.^-k - rn.^k).*rq./ ... 
            (rq.^2 - rn.^2).*(rn.^k./rq.^k - rq.^k./rn.^k);
    rnc = rp./rn.*(rp.^-k - rp.^k)./(rn.^-k - rn.^k);
    rqc_rnzero = rp.*rq.^(k-1).*(rp.^-k - rp.^k)./k; % rn = 0 treated differently!
    rqc(rn==0) = rqc_rnzero(rn==0);
    rqc(isnan(rqc)) = 0;    % remove NaN
    rnc(isnan(rnc)) = 0;
    % oxygen perturbation effect
    oxypert = -rp.^(-k-1).*rq.^(k+1)*mu_prol.*rqc - ...  
               rp.^(-k-1).*rn.^(k+1)*mu_die.*rnc;
    % get Lambda(k)
    Lambda(k, :) = mu_prol*0.5*(mu_die/mu_prol*rn.^2.*(1+k) + rq.^2.*(1+k) ...
                         + rp.^2.*(1-k))./rp.^2 - sigma.*k.*(k.^2-1)./rp.^3 ...
                         + oxypert;
    sigma_stable(k, :) = (1+oxypert)./(k.*(k.^2-1)./rp.^3); % at stationary state
end
end