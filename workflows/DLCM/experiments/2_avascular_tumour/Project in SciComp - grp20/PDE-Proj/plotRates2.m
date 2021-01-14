%% Plot the rates

% Plot the rate_sdof as bars
bartender = zeros(1,size(inspect_rates,2));
for kk = 1:length(bartender)
    rate_sdof = inspect_rates{1,kk};
    sdof_m_ = inspect_rates{2,kk};
    idof_ = inspect_rates{3,kk};
    bartender(kk) = sum(rate_sdof(idof_));
end
b = bar(bartender,'LineStyle','none');