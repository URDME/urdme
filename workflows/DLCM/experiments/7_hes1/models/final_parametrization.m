%FINAL_PARAMETRIZATION Final parameterisation of parameters.
%   Check parameterisation of perturbed concentrations and mu
%   parameters followed by finding perturbed alpha parameters and
%   fitting to lognormal distributions.

% G. Menz 2024-10-22

%% (0) Set up

% find 68% confidence interval for lognormal distribution, so seek
% parameters x([1 2]) s.t.
fun = @(x)logninv([(1-0.68)/2 1-(1-0.68)/2],x(1),x(2));

NMC = 1000;

%% (1) Check distribution of perturbed concentrations

% assumption: lognormally distributed with 5% noise
conc = hes1_conc;
conc_pert = hes1_conc(NMC);

% start from a given relative perturbation level, i.e. 5% for concentrations
xx = [log(conc) repmat(0.05,numel(conc),1)];

% find 68% confidence intervals using inverse lognormal distribution
conc_cis = zeros(numel(conc),2);
for i = 1:numel(conc)
  conc_cis(i,:) = fun(xx(i,:));
end

% plot histograms, means and CIs of all perturbed concentrations
figure(1), clf,
t = tiledlayout(1,5);
title(t,'Perturbed concentrations from hes1_conc','Interpreter','none')

nexttile;
histogram(conc_pert(1,:),'Normalization','pdf')
hold on
xline(conc(1),'r','LineWidth',2)
hold on
xline(conc_cis(1,:),'b','LineWidth',2)
xlabel('Dll1 conc')

nexttile;
histogram(conc_pert(2,:),'Normalization','pdf')
hold on
xline(conc(2),'r','LineWidth',2)
hold on
xline(conc_cis(2,:),'b','LineWidth',2)
xlabel('Notch conc')

nexttile;
histogram(conc_pert(3,:),'Normalization','pdf')
hold on
xline(conc(3),'r','LineWidth',2)
hold on
xline(conc_cis(3,:),'b','LineWidth',2)
xlabel('Hes1 mRNA conc')

nexttile;
histogram(conc_pert(4,:),'Normalization','pdf')
hold on
xline(conc(4),'r','LineWidth',2)
hold on
xline(conc_cis(4,:),'b','LineWidth',2)
xlabel('Hes1 protein conc')

nexttile;
histogram(conc_pert(5,:),'Normalization','pdf')
hold on
xline(conc(5),'r','LineWidth',2)
hold on
xline(conc_cis(5,:),'b','LineWidth',2)
xlabel('Ngn2 conc')

%% (2) Check distribution of perturbed mu parameters

% find parameters (both perturbed and unperturbed)
par = hes1_params;
par_pert = hes1_params(NMC);

% for muD and muN we assume a lognormal distribution with 10% noise so we 
% find 68% CIs similarly to above
xx1 = [log(par.muD) 0.1; log(par.muN) 0.1];
mu_cis = zeros(5,2);
mu_cis(1,:) = fun(xx1(1,:));
mu_cis(2,:) = fun(xx1(2,:));

% for muM, muP and mun we are given 68% CIs in references
mu_cis(3,:) = [log(2)/25.8 log(2)/22.4]; % half-life: 24.1 pm 1.7 min
mu_cis(4,:) = [log(2)/25.4 log(2)/19.2]; % half-life: 22.3 pm 3.1 min
mu_cis(5,:) = [log(2)/24.1 log(2)/19.7]; % half-life: 21.9 pm 2.2 min

% find distributions from these CIs
mu_mus = zeros(5,1);
mu_sigmas = zeros(5,1);
for i = 1:5
  % starting guess: middle of interval and its width:
  M = mean(mu_cis(i,:));
  V = diff(mu_cis(i,:))^2;
  % transfer (M,V) to underlying normal distribution (see lognrnd)
  MU = log(M^2 / sqrt(V+M^2));
  SIGMA = sqrt(log(V/M^2 + 1));

  % solve
  xx = fsolve(@(x)fun(x)-mu_cis(i,:),[MU SIGMA],optimset('Display','off'));
  mu_mus(i) = xx(1);
  mu_sigmas(i) = xx(2);
end

% plot histograms, pdfs, means and CIs for all mu parameters
figure(2), clf,
t = tiledlayout(1,5);
title(t,'Perturbed mus from hes1_params','Interpreter','none')

nexttile;
histogram(par_pert.muD,'Normalization','pdf')
hold on
xline(par.muD,'r','LineWidth',2)
hold on
xline(mu_cis(1,:),'b','LineWidth',2)
hold on
x1 = linspace(0.05,0.1);
p1 = lognpdf(x1,mu_mus(1),mu_sigmas(1));
plot(x1,p1,'k','LineWidth',2)
xlabel('$\mu_D$','Interpreter','latex')

nexttile;
histogram(par_pert.muN,'Normalization','pdf')
hold on
xline(par.muN,'r','LineWidth',2)
hold on
xline(mu_cis(2,:),'b','LineWidth',2)
hold on
x2 = linspace(0.06,0.12);
p2 = lognpdf(x2,mu_mus(2),mu_sigmas(2));
plot(x2,p2,'k','LineWidth',2)
xlabel('$\mu_N$','Interpreter','latex')

nexttile;
histogram(par_pert.muM,'Normalization','pdf')
hold on
xline(par.muM,'r','LineWidth',2)
hold on
xline(mu_cis(3,:),'b','LineWidth',2)
hold on
x3 = linspace(0.02,0.04);
p3 = lognpdf(x3,mu_mus(3),mu_sigmas(3));
plot(x3,p3,'k','LineWidth',2)
xlabel('$\mu_M$','Interpreter','latex')

nexttile;
histogram(par_pert.muP,'Normalization','pdf')
hold on
xline(par.muP,'r','LineWidth',2)
hold on
xline(mu_cis(4,:),'b','LineWidth',2)
hold on
x4 = linspace(0.02,0.048);
p4 = lognpdf(x4,mu_mus(4),mu_sigmas(4));
plot(x4,p4,'k','LineWidth',2)
xlabel('$\mu_P$','Interpreter','latex')

nexttile;
histogram(par_pert.mun,'Normalization','pdf')
hold on
xline(par.mun,'r','LineWidth',2)
hold on
xline(mu_cis(5,:),'b','LineWidth',2)
hold on
x5 = linspace(0.02,0.042);
p5 = lognpdf(x5,mu_mus(5),mu_sigmas(5));
plot(x5,p5,'k','LineWidth',2)
xlabel('$\mu_n$','Interpreter','latex')

%% (3) Parameterise alpha parameters using perturbed concentrations and mus

% load data or use parametrization.m for parameterisation
if ~exist('../data/alpha_parameterization.mat','file')
  final = true;
  parametrization; % takes a while to run for high NMC
else
  load('../data/alpha_parameterization.mat')
end

% initialise vectors
dists = cell(5,1);
alpha_mus = zeros(5,1);
alpha_sigmas = zeros(5,1);
means = zeros(5,1);
alpha_cis = zeros(5,2);

% fit lognormal distribution for all alphas
for i = 1:5
  dists{i} = fitdist(alphas(i,:)','Lognormal');
  pv = dists{i}.ParameterValues;
  alpha_mus(i) = pv(1);
  alpha_sigmas(i) = pv(2);

  % find mean and CI
  means(i) = mean(dists{i});
  alpha_cis(i,:) = icdf(dists{i},[(1-0.68)/2 1-(1-0.68)/2]);
end

%% (4) Check alpha parameterisation and the found distributions
figure(3),clf,
t = tiledlayout(1,5);
title(t,'Final alphas after parameterisation')

% plot histogram, distribution, mean and CI for all alphas
nexttile;
histogram(alphas(1,:),20,'Normalization','pdf')
hold on
x1 = linspace(0.01,0.03);
p1 = lognpdf(x1,dists{1}.mu,dists{1}.sigma);
plot(x1,p1,'k','LineWidth',2); hold on,
xline(means(1),'r','LineWidth',2); hold on,
xline(alpha_cis(1,:),'b','LineWidth',2); 
xlabel('$\alpha_D$','Interpreter','latex')

nexttile;
histogram(alphas(2,:),20,'Normalization','pdf')
hold on
x2 = linspace(3.8,9.2);
p2 = lognpdf(x2,dists{2}.mu,dists{2}.sigma);
plot(x2,p2,'k','LineWidth',2); hold on,
xline(means(2),'r','LineWidth',2); hold on,
xline(alpha_cis(2,:),'b','LineWidth',2);
xlabel('$\alpha_N$','Interpreter','latex')

nexttile;
histogram(alphas(3,:),20,'Normalization','pdf')
hold on
x3 = linspace(0.01,0.03);
p3 = lognpdf(x3,dists{3}.mu,dists{3}.sigma);
plot(x3,p3,'k','LineWidth',2); hold on,
xline(means(3),'r','LineWidth',2); hold on,
xline(alpha_cis(3,:),'b','LineWidth',2);
xlabel('$\alpha_M$','Interpreter','latex')

nexttile;
histogram(alphas(4,:),20,'Normalization','pdf')
hold on
x4 = linspace(0.05,0.24);
p4 = lognpdf(x4,dists{4}.mu,dists{4}.sigma);
plot(x4,p4,'k','LineWidth',2); hold on,
xline(means(4),'r','LineWidth',2); hold on,
xline(alpha_cis(4,:),'b','LineWidth',2); 
xlabel('$\alpha_P$','Interpreter','latex')

nexttile;
histogram(alphas(5,:),20,'Normalization','pdf')
hold on
x5 = linspace(3.2e-3,7.2e-3);
p5 = lognpdf(x5,dists{5}.mu,dists{5}.sigma);
plot(x5,p5,'k','LineWidth',2); hold on,
xline(means(5),'r','LineWidth',2); hold on,
xline(alpha_cis(5,:),'b','LineWidth',2);
xlabel('$\alpha_n$','Interpreter','latex')

%% (5) Check alpha parameterisation implemented in hes1_params

% plot histogram, distribution, mean and CI for all alphas
figure(4), clf,
t = tiledlayout(1,5);
title(t,'Perturbed alphas from hes1_params','Interpreter','none')

nexttile;
histogram(par_pert.alphaD,'Normalization','pdf')
hold on
xline(par.alphaD,'r','LineWidth',2)
hold on
xline(alpha_cis(1,:),'b','LineWidth',2)
hold on
plot(x1,p1,'k','LineWidth',2);
xlabel('$\alpha_D$','Interpreter','latex')

nexttile;
histogram(par_pert.alphaN,'Normalization','pdf')
hold on
xline(par.alphaN,'r','LineWidth',2)
hold on
xline(alpha_cis(2,:),'b','LineWidth',2)
hold on
plot(x2,p2,'k','LineWidth',2);
xlabel('$\alpha_N$','Interpreter','latex')

nexttile;
histogram(par_pert.alphaM,'Normalization','pdf')
hold on
xline(par.alphaM,'r','LineWidth',2)
hold on
xline(alpha_cis(3,:),'b','LineWidth',2)
hold on
plot(x3,p3,'k','LineWidth',2);
xlabel('$\alpha_M$','Interpreter','latex')

nexttile;
histogram(par_pert.alphaP,'Normalization','pdf')
hold on
xline(par.alphaP,'r','LineWidth',2)
hold on
xline(alpha_cis(4,:),'b','LineWidth',2)
hold on
plot(x4,p4,'k','LineWidth',2);
xlabel('$\alpha_P$','Interpreter','latex')

nexttile;
histogram(par_pert.alphan,'Normalization','pdf')
hold on
xline(par.alphan,'r','LineWidth',2)
hold on
xline(alpha_cis(5,:),'b','LineWidth',2)
hold on
plot(x5,p5,'k','LineWidth',2);
xlabel('$\alpha_n$','Interpreter','latex')
