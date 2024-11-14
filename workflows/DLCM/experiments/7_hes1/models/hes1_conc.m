function conc = hes1_conc(NMC)
%HES1_CONC Concentrations of all constitutents in the Hes1-model.
%   The values are determined as described in the associated
%   manuscript (see README).
% 
%   Examples:
%     % mean concentrations
%     conc = hes1_conc;
%
%     % 50 perturbed concentrations
%     conc_perturbed = hes1_conc(50);
%
%   See also HES1_PARAMS, HES1RED_PARAMS.

% G. Menz 2024-09-17

% units: uM
conc = zeros(5,1);
conc(1) = 0.0135;   % "D" Dll1 conc
conc(2) = 0.925;    % "N" Notch conc
conc(3) = 0.0613;   % "M" Hes1 mRNA conc
conc(4) = 0.269;    % "P" Hes1 protein conc
conc(5) = 0.0505;   % "n" Ngn2 conc

% optionally generate NMC perturbed concentrations
if nargin > 0
  % ad hoc assumption: lognormal distribution with ~5% noise
  conc =  conc.*exp(0.05*randn(5,NMC));
end
