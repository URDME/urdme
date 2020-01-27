function current = getCurrent(compartment,t)
%GETCURRENT Current in compartment at a given time.
%   CURRENT = GETCURRENT(COMPARTMENT,T) is designed to be called by
%   COMSOL and to return the current of the compartment in question.
%
%   To plot it in COMSOL, create a line-plot in the results of COMSOL
%   and specify compartment and use time on the x-axis.

% A. Senek 2017-07-12

% interpolation data, to be computed offline
persistent I;
if isempty(I)
  load('I.mat');
end
% 1st column is time in ms, next follows the compartments

% current is in nA, hence 1e-9
current = 1e-9*interp1q(I(:,1),I(:,compartment),1000*t');
% in case the above is no longer supported:
%current = 1e-9*interp1(I(:,1),I(:,compartment),1000*t');
