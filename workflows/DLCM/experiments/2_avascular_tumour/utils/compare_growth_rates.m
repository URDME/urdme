% COMPARE_GROWTH_RATES
% This code uses experimental data from AvascularTumor_updated.m to
% estimate amplitude growth of modes of a perturbed tumor boundary.
% The analytical value for the amplitude growth rate based on the
% dispersion relation obtained from linear stability analysis and the tumor
% model region volume characteristics: proliferating, quiescent and, 
% necrotic region.
% 
% The estimated experimental growth is based on MSE fit, assuming an
% exponentially growing/decaying perturbation mode.
%
% Experimental data from 'growth_rate_comparison_exp.mat'.
%
% E. Blom 2024-02-01

% Run 1st section to evaluate Lambda from the data
% Run 2nd section to directly get the Figure from the article


%% Evaluate and estimate Lambda from data
% exponential growth function
modelfun = @(b, t)b(1)*exp(b(2)*t);

% load data from all experiments
filename = "growth_rate_comparison_exp.mat";
load(filename)  % load cell-struct A with all data
% load parameters from 1st experiment
lambda = cell_data{1,1}.all_parameters(1);
kappa_prol = cell_data{1,1}.all_parameters(2);
mu_prol = cell_data{1,1}.all_parameters(3);
kappa_death = cell_data{1,1}.all_parameters(4);
mu_death = cell_data{1,1}.all_parameters(5);
mu_degrade = cell_data{1,1}.all_parameters(6);
sigma = cell_data{1,1}.all_parameters(7);

% time vector tspan assumed same for all experiments
tspan = cell_data{1,1}.tspan;

% estimate amplitude of boundary perturbation modes 1:J
[J, Ns] = size(cell_data);         % k = 1,...,J, Ns = number of samples per k

T = numel(cell_data{1,1}.tspan);   % nr time steps in simulation data
amp_k = zeros(T-1, J, Ns); % observed perturbation mode amplitude


for k = 1:J % for all modes
for n = 1:Ns % for all experiments
    
    % load experiment data
    BFTsave = cell_data{k,n}.BFTsave;
    rp = cell_data{k,n}.rp;
    rq = cell_data{k,n}.rq;
    rn = cell_data{k,n}.rn;
    % find mode amplitude for mode k at each time step
    for i = 2:numel(tspan)
        amp_k(i-1,k, n) = BFTsave{i}(2,k+1);
    end

    % find analytical growth rates
    args = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
    'rp',rp,'rq',rq, 'rn',rn,'sigma',sigma,'modes',k) ;
    [Lambda, sigma_req] = get_growth_factors(args);
    Lambda = Lambda./2.7;               % map to DLCM model units
    Lambda_k(:, k, n) = Lambda(k,:);    % save all analytical values

end
% find growth rates from SAMPLE means
Lambda_k_mean = mean(Lambda_k(:,k,:), 3);   % analytical
amp_k_mean = mean(amp_k(:,k,:), 3);         % experimental
figure(2); clf;
plot(Lambda_k_mean);                        % plot SAMPLE means...
xlabel('Time', 'Interpreter','latex')
ylabel('$\Lambda(k)$', 'Interpreter','latex')
figure(1); clf;
plot(tspan(2:end), amp_k_mean(:));
xlabel('Time', 'Interpreter','latex')
ylabel('Observed mode amplitude', 'Interpreter','latex')
% ...and ask user to identify approximately steady time intervals
ta = input('Starting point for interval of ~constant Lambda: ');
tb = input('End point for interval of ~constant Lambda: ');

% analytical mean during time interval
LAM = mean(Lambda_k_mean(ta:tb)); % after 'burn-in' period
LAM_std = std(Lambda_k_mean(ta:tb));

% experimental mean during time interval:
% put all data in one vector...
T = tspan(ta+1:tb+1)'.*ones(1,20);
Y = log(amp_k(ta:tb,k,:));
T = T(:);
Y = Y(:);

% ... and fit log of observed amplitude for observed growth rate and
% uncertainty estimate
pfitdata = fit(T(:),Y(:),'poly1')
% and show analytical value with std
disp("analytical Lambda(" + k + "), mean: " + LAM + ", std: " + LAM_std)
end

return;

%% Plot
% Experimental and analytical Lambda estimated from the data from 8*20 runs
% using Section 1 above
lam_est = [0.074, 0.17, 0.16, 0.119, 0.0374, -0.139, -0.37, -0.42];
lam_est_std = [0.0018, 0.0017, 0.002, 0.0027, 0.009, 0.019, 0.016, 0.019];
LAM = [0.026, 0.14, 0.158, 0.126, 0.063, -0.072, -0.227, -0.43];
LAM_std = [0.0006, 0.002, 0.003, 0.004, 0.004, 0.008, 0.023, 0.024];

% plot
figure(2), clf;
k = 1:8;
errorshade(k, lam_est - lam_est_std, lam_est + lam_est_std)
hold on; plot(k, lam_est, '^')
errorbar(k, LAM, LAM_std, '-.')
xlabel('Mode k', 'Interpreter','latex')
ylabel('$\Lambda(k)$ $[\mu_{\mathrm{prol}}]$', 'Interpreter','latex')
axis([k(1)-0.3,k(end)+0.3, min(LAM)-0.05, max(LAM)+0.05])
grid on
legend('experimental', 'analytical')
set(gcf,'Position',[100 100 340/(sqrt(2)) 240/(sqrt(2))]);
set(gcf,'PaperPositionMode','auto');
set(gca, 'fontname', 'Roman', 'FontSize', 9.0)
%$$$ uncomment to save:
%$$$ exportgraphics(gca,'Fig_new.pdf')