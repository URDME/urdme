%% Plot the rates

% Define the different colors for the rates to be the same always
keySet = {'moveb','moveb1','moveb2','moves','birth','death','degrade'};
valueSet = {'#0072BD','#4DBEEE','#A2142F','#D95319','#EDB120','#7E2F8E','#77AC30'};
RateColors = containers.Map(keySet,valueSet);

% Get the name tag of each name (used in the legend)
rate_names = fieldnames(Ne);

% Plot the rates as stacked bars (non- or normalized)
rates = abs(inspect_rates);
% rates = rates./sum(rates,1); % normalizing
b = bar(rates','stacked','LineStyle','none');

% Apply the appropriate color
for kk = 1:length(rate_names)
    b(kk).FaceColor = RateColors(rate_names{kk});
end

% Plot settings
grid on;
title(sprintf('Relative and normalized rates \n alpha = %d', alpha));
xlabel('time')
ylabel('rates')
% ticks = 
set(gca, 'XTick', linspace(1,length(tspan),7))
set(gca, 'XTickLabel', round(linspace(1,tspan(end),7)))
% ylim([0 1.5]);
legend(rate_names);