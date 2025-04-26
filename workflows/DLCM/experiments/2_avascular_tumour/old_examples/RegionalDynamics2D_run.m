% runs RegionalDynamics2D and displays result
% includes calculations of perturbation factors during growth

% dependencies: 
%   - RegionalDynamics2D.m
%   - get_growth_factors.m
%   

%% Numerical solution, radial growth dynamics
% model parameters and initial condition
kappa_prol = 0.94;
kappa_death = 0.93;
mu_prol = 1;
mu_death = 1.35;
lambda = 1.15;
c_out = 1;
init_radius = 0.1;
init_volume = init_radius^2*pi;
Tend = 12;
sigma = 0e-6;      % surface tension parameter

% simulate 1D dynamics
par = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
             'kappa_prol',kappa_prol,'kappa_death',kappa_death, ...
             'lambda',lambda,'c_out',c_out);
reg_char = RegionalDynamics2D(par, init_volume, linspace(0,Tend,1000));                            
Vp = reg_char.Vp;
Vq = reg_char.Vq;
Vn = reg_char.Vn;
tspan = reg_char.tspan;
radial_discriminant = reg_char.discr;

% confirm that regions satisy invariants K_prol, K_death
K_prol = 4*(1-kappa_prol)/lambda;
K_death = 4*(1-kappa_death)/lambda;
figure(3), clf,
rq = sqrt(Vq/pi); rp = sqrt(Vp/pi); rn = sqrt(Vn/pi);
semilogy(tspan, ...
         abs(-(rp.^2.*log(rp.^2) - rn.^2.*log(rq.^2) + (rq.^2 - rp.^2))-K_prol));
hold on,
semilogy(tspan, ...
         abs(-(rp.^2.*log(rp.^2) - rn.^2.*log(rn.^2) + (rn.^2 - rp.^2))-K_death));
legend('Error K_{prol}','Error K_{death}');
xlabel('time')

% time intervals for plotting markers
interval = 10;
tint = round(length(tspan)/interval);
t_idx = sort([2:tint:numel(tspan)]);

% plot the regional characteristics
figure(1), clf
p1 = plot(tspan(t_idx), Vp(t_idx), '*', 'LineWidth', 1);
hold on
p2 = plot(tspan(t_idx), Vq(t_idx), '.', 'Marker', '.', 'MarkerSize', 15);
plot(tspan(t_idx), Vn(t_idx), 'kx', 'LineWidth', 1)
p1.Color = graphics_color('vermillion');
p2.Color = graphics_color('bluish green');
p3 = plot(tspan, Vp, '-', 'LineWidth', 1);
p4 = plot(tspan, Vq, '-', 'LineWidth', 1);
plot(tspan, Vn, 'k-', 'LineWidth', 1);
p3.Color = graphics_color('vermillion');
p4.Color = graphics_color('bluish green');
%plot(tspan(1:end-1), discriminant, '--', 'LineWidth', 2)
%plot(tspan(1:end-1), sonr + 0.*tspan(1:end-1), 'k--')
xlabel('Time [$\mu_{prol.}^{-1}$]', 'Interpreter','latex')
ylabel('Volume [$V_R$]', 'Interpreter','latex')
axis( [0,tspan(end), 0,0.4])

%yyaxis right
%plot(tspan(t_idx), radial_discriminant(t_idx), '--', 'LineWidth', 2)
%ylabel('$\Delta_{\theta} [\mu_{prol.}]$', 'Interpreter', 'Latex')

legend('$V_p$', '$V_q$', '$V_n$', 'Interpreter', 'Latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340/(sqrt(2)) 240/(sqrt(2))]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
% uncomment to save
% exportgraphics(gca,'regionaldynamics_ex.pdf')

% plot analytical perturbation coefficients
modes = [1 3, 9]; ...        % for these modes
args = struct('mu_prol',mu_prol,'mu_death',mu_death, ...
             'rp',sqrt(Vp/pi),'rq',sqrt(Vq/pi), ...
             'rn',sqrt(Vn/pi),'sigma',sigma,'modes',modes) ;
Lambda = get_growth_factors(args);
figure(2), clf,
% plot markers
plot(tspan(t_idx), Lambda(modes(1),t_idx), 'o', 'LineWidth', 1)
hold on
plot(tspan(t_idx), Lambda(modes(2),t_idx), '^', 'LineWidth', 1)
plot(tspan(t_idx), Lambda(modes(3),t_idx), 'square', 'LineWidth', 1)
%plot(tspan(2:tint:end-1), Lambda(modes(4),2:tint:end), '-pentagram', 'LineWidth', 1, 'MarkerSize', 10)
%hold on; 
% plot positive term for k=1 (markers)
 plot(tspan(t_idx), mu_prol*0.5*(mu_death/mu_prol*Vn(t_idx).*(1+1) + Vq(t_idx).*(1+1) ...
                          + Vp(t_idx).*(1-1))./Vp(t_idx), 'k+', 'MarkerSize', 10)
 % % plot negative terms for k=1 
 % yyaxis right
 plot(tspan(t_idx), -(Lambda(1,t_idx)' - mu_prol*0.5*(mu_death/mu_prol*Vn(t_idx).*(1+1) + Vq(t_idx).*(1+1) ...
                          + Vp(t_idx).*(1-1))./Vp(t_idx)), 'k_', 'MarkerSize', 10)
 % plot full lines
 plot(tspan, Lambda(modes(1),:), '-', 'LineWidth', 1,  'color', [0 0.4470 0.7410])
 hold on
 plot(tspan, Lambda(modes(2),:), '-', 'LineWidth', 1,  'color', [0.8500 0.3250 0.0980])
 plot(tspan, Lambda(modes(3),:), '-', 'LineWidth', 1,  'color', [0.9290 0.6940 0.1250])
 plot(tspan, mu_prol*0.5*(mu_death/mu_prol*Vn.*(1+1) + Vq.*(1+1) ...
                          + Vp.*(1-1))./Vp, 'k:')
 plot(tspan, -(Lambda(1,:)' - mu_prol*0.5*(mu_death/mu_prol*Vn.*(1+1) + Vq.*(1+1) ...
                          + Vp.*(1-1))./Vp), 'k:')

 xlabel('Time [$\mu_{prol.}^{-1}$]', 'Interpreter','latex')
 ylabel('$\Lambda(k)$', 'Interpreter','latex')
% 
 axis([0,Tend, 0, 1.1])
 hold on
 
% axis([0,tspan(end), -1 0])
% ylabel('negative contribution, $k=1$', 'Interpreter','latex')
legend('$k=1$', '$k=3$', '$k=9$', 'Interpreter','latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[400 100 340/(sqrt(2)) 240/(sqrt(2))]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
% uncomment to save
% exportgraphics(gca,'lambdadynamics_ex.pdf')
