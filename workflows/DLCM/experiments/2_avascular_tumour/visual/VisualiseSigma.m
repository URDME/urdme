% Plots Lambda(k) versus k, given tumor regional sizes Vp, Vq, Vn

% E. Blom 2023-02-01

modes = 1:15;
sigma_range = [0, 5*1e-5, 5*1e-4, 5*1e-3];

% Tumor geometric state
Vp = 0.31;  % tumor size
Vq = 0.24;  % quiescent-proliferative interface enclosed region
Vn = 0.052;   % necrotic-quiescent interface enclosed region
 
args = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
             'rp',sqrt(Vp/pi),'rq',sqrt(Vq/pi), ...
             'rn',sqrt(Vn/pi),'sigma',sigma_range,'modes',modes) ;
Lambda = get_growth_factors(args);

plot(modes, Lambda(modes,1), '-o', 'Linewidth', 1)
hold on
plot(modes, Lambda(modes,2), '-^', 'Linewidth', 1)
plot(modes, Lambda(modes,3), '-square', 'Linewidth', 1)
plot(modes, Lambda(modes,4), '-pentagram', 'Linewidth', 1, 'MarkerSize', 10)
plot(modes, 0.*modes, 'k--')
%lgd = legend(sprintf("%0.0E", sigma(1)),...
%    sprintf("%0.0E", sigma(2)), ...
%    sprintf("%0.0E", sigma(3)), ...
%    sprintf("%0.0E", sigma(4)), 'Interpreter', 'Latex');
legend('$\sigma = 0$',...
'$5\cdot 10^{-5}$', ...
'$5\cdot 10^{-4}$', ...
'$5\cdot 10^{-3}$', 'Interpreter', 'Latex');
warning('legend set manually and values might not match with input')

xlabel('Mode k', 'Interpreter','latex')
ylabel('$\Lambda(k)$', 'Interpreter','latex')
axis([modes(1)-1,modes(end)+1, -0.3, 1.1])
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340/(sqrt(2)) 240/(sqrt(2))]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
% uncomment to save
%exportgraphics(gca,'sigmadynamics_ex1.pdf')