function params = hes1_params(NMC)
%HES1_PARAMS Parameter values for the Hes1@5 model.
%
%   Examples:
%     % call with no input to get original parameter values
%     par = hes1_params;
%
%     % get 50 perturbed parameters according to the manuscript (see
%     % README):
%     pars = hes1_params(50);
%
%   See also HES1_CONC, HES1RED_PARAMS.

% G. Menz 2024-02-26

% parameter values scaled to fit wanted concentrations (see
% parametrization for more information)
params.alphaD = 0.018499886514315;
params.alphaN = 5.983559710192386;
params.alphaM = 0.017418330115810;
params.alphaP = 0.139114902363266;
params.alphan = 0.004853715617249;

% factor by which muD and muN are multiplied: if factor = 5, we have
% 80% efficiency, i.e., 80% is used, while 20% is degraded
params.muD = log(2)/50*5;
params.muN = log(2)/40*5;
params.muM = log(2)/24.1;
params.muP = log(2)/22.3;
params.mun = log(2)/21.9;

% Hill function parameters
params.k = 1;
params.h = 4;
params.KM = 0.05; 
params.Kn = 0.03;

% generate perturbed parameters
if nargin > 0
  % alphas fitted to lognormal distributions in final_parametrization
  alphas_mus = [-3.997551210995813; ...
                1.782188522311404; ...
                -4.056082749832997; ...
                -1.985690530006482; ...
                -5.333996532887211];
  alphas_sigmas = [0.122967717046486; ...
                   0.116851511000996; ...
                   0.108171877201908; ...
                   0.162698975945510; ...
                   0.109414552906170];
  params.alphaD = lognrnd(alphas_mus(1),alphas_sigmas(1), ...
                          1,NMC);
  params.alphaN = lognrnd(alphas_mus(2),alphas_sigmas(2), ...
                          1,NMC);
  params.alphaM = lognrnd(alphas_mus(3),alphas_sigmas(3), ...
                          1,NMC);
  params.alphaP = lognrnd(alphas_mus(4),alphas_sigmas(4), ...
                          1,NMC);
  params.alphan = lognrnd(alphas_mus(5),alphas_sigmas(5), ...
                          1,NMC);

  % ad hoc 10% noise for muN and muD
  params.muD = params.muD*exp(0.1*randn(1,NMC));
  params.muN = params.muN*exp(0.1*randn(1,NMC));

  % assumption: lognormal distribution of muM with std 1.7 fitting distn
  % to 68% CI gives the following parameters
  mu_muM = -3.546227554841154;
  sigma_muM = 0.071049163683391;
  params.muM = lognrnd(mu_muM,sigma_muM,1,NMC);

  % assumption: lognormal distribution of muP with std 3.1 fittting
  % distn to 68% CI gives the following parameters
  mu_muP = -3.461291669344817;
  sigma_muP = 0.140661642709813;
  params.muP = lognrnd(mu_muP,sigma_muP,1,NMC);

  % assumption: lognormal distribution of mun with S.E.M. 2.2 fitting
  % distn to 68% CI gives the following parameters
  mu_mun = -3.447927657936440;
  sigma_mun = 0.101350529532195;
  params.mun = lognrnd(mu_mun,sigma_mun,1,NMC);

  % Hill-function assumed given
  params.k = repmat(params.k,1,NMC);
  params.h = repmat(params.h,1,NMC);
  params.KM = repmat(params.KM,1,NMC);
  params.Kn = repmat(params.Kn,1,NMC);
end
