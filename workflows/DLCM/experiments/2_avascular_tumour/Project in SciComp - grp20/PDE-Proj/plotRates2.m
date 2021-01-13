%% Plot the rates

% Plot the rate_sdof as bars
figure;
hold on;
bartender = zeros(1,size(inspect_rates,2));
for kk = 1:length(bartender)
    rate_sdof = inspect_rates{1,kk};
    sdof_m_ = inspect_rates{2,kk};
    idof_ = inspect_rates{3,kk};
    bartender(kk) = sum(rate_sdof(idof_));
end
b = bar(bartender,'LineStyle','none');
% Plot settings

% title(sprintf('Relative and normalized rates \n alpha = %d', alpha));
% xlabel('time')
% ylabel('rates')
% % ticks = 
% set(gca, 'XTick', linspace(1,length(tspan),7))
% set(gca, 'XTickLabel', round(linspace(1,tspan(end),7)))
% % ylim([0 1.5]);
% legend(rate_names);