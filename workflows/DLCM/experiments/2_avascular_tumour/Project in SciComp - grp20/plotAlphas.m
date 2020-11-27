% Plot alphas

% Load files
U1 = load('saveData_alpha1e-01_2020-11-24-9-49-14.mat','U').U;
U4 = load('saveData_alpha4e-01_2020-11-24-9-54-18.mat','U').U;
U8 = load('saveData_alpha8e-01_2020-11-24-10-10-43.mat','U').U;

alphas = [0.1 0.2 0.4 0.6 0.8 1];
N_cells = [sum(abs(U1)) 1850 sum(abs(U4)) 2019 sum(abs(U8)) 2321];

IC = 5;
Tend = 100;

% Plot number of cells / alpha
plot(alphas, N_cells, '-*', 'Linewidth', 1.5);
xlabel('\alpha');
ylabel('Number of cells');
title(sprintf('IC = %d, Time = %d', ...
        IC, Tend));
set(gca,'FontSize',14);