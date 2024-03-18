% Compare simulated 2D tumor geometric characteristics
% with corresponding analytical radial tumor growth characteristics.
%
% *** Run after running AvascularTumorPDE.m ***

% E. Blom 2022-12-01

% Dependencies:
%   - RegionalDynamics2D

%% Regional volume characteristics

plotleft = false; 
show_full_time = false; % show until final tspan or tstamp

% filtering roundness only
sgf_wind = 21;   % savitzky-golay filter window
sgf_deg = 3;

c_out = 1;
mu_death_eff = 1.35; %1.0
lambda_eff = 1.15;   %1.1
% find approximated 1D analytical solution
init_volume = rp(2)^2*pi;
par = struct('mu_prol',mu_prol,'mu_death',mu_death_eff, ...
             'kappa_prol',kappa_prol,'kappa_death',kappa_death, ...
             'lambda',lambda_eff,'c_out',c_out);
reg_char = RegionalDynamics2D(par, init_volume, linspace(0,Tend,1000));                            

Vp = reg_char.Vp;
Vq = reg_char.Vq;
Vn = reg_char.Vn;
radial_discriminant = reg_char.discr;
tspan_analytical = reg_char.tspan;

% plot in same units of time
tscale = 2.7;       % timescale factor
tspan_analytical_scaled = tscale.*tspan_analytical;

if PDE
    % rescale the simulated units of time
    tspan_scaled = tscale.*tspan;
else
    tspan_scaled = tspan;
end

figure(4)
colororder({'k', '[0 0.4470 0.7410]' })

% Plot 2D numerical dynamics
% % !!! JUST FOR FIG 4.1
% rp = sqrt(Vp/pi);
% rq = sqrt(Vq/pi);
% rn = sqrt(Vn/pi);
% % !!!
tint = round(length(tspan_scaled)/20); 
p5 = plot(tspan_scaled([2:tint:end end]), rp([2:tint:end end]).*rp([2:tint:end end])*pi, '--*', 'LineWidth', 1);
hold on
p6 = plot(tspan_scaled([2:tint:end end]), rq([2:tint:end end]).*rq([2:tint:end end])*pi, '--.', 'LineWidth', 1, 'MarkerSize', 15);
plot(tspan_scaled([2:tint:end end]), rn([2:tint:end end]).*rn([2:tint:end end])*pi, 'k--x', 'LineWidth', 1);
p5.Color = graphics_color('vermillion');
p6.Color = graphics_color('bluish green');

% Plot approximated 1D analytical dynamics
p1 = plot(tspan_analytical_scaled, Vp, 'LineWidth', 1);
p2 = plot(tspan_analytical_scaled, Vq, 'LineWidth', 1);
plot(tspan_analytical_scaled, Vn, 'k', 'LineWidth', 1)
%p3 = plot(tspan_analytical(1:end-1), discriminant, '--', 'LineWidth', 2);
%p4 = plot(tspan_analytical(1:end-1), sonr + 0.*tspan_analytical(1:end-1));
p1.Color = graphics_color('vermillion');
p2.Color = graphics_color('bluish green');
%p3.Color = graphics_color('orange');
%p4.Color = graphics_color('blue');

% selected frames for spatial visualisation
t_frames = [tspan_scaled(tstamps(1)), tspan_scaled(tstamps(2)), tspan_scaled(tstamps(3))]; 
plot([t_frames(1) t_frames(1)], [min(Vn) max(Vp)*1.8], 'k-.', 'LineWidth', 1.25)
plot([t_frames(2) t_frames(2)], [min(Vn) max(Vp)*1.8], 'k-.', 'LineWidth', 1.25)
plot([t_frames(3) t_frames(3)], [min(Vn) max(Vp)*1.8], 'k-.', 'LineWidth', 1.25)

%title('Tumor regional characteristics',  'Interpreter','latex')

if plotleft
    xlabel('Time [$\mu_{\mathrm{prol}}^{-1}$]', 'Interpreter','latex')
    ylabel('Volume [$\pi R^2$]', 'Interpreter','latex')
else
    xlabel('Time [$\mu_{\mathrm{prol}}^{-1}$]', 'Interpreter','latex')
    set(gca,'ytick',[])
end

if show_full_time
    axis([0, tspan_scaled(end), 0, 0.5])
else % show until final tstamp
    axis([0, tspan_scaled(tstamps(end)), 0, 0.5])
end

yyaxis right
% filter roundness
% roundness not defined after tumor splits
if plotleft
    final_idx = find(diff(roundness) < -10);  % different noise levels
else
    final_idx = find(diff(roundness) < -1);   % different noise levels
end
if isempty(final_idx)
    final_idx = numel(tspan_scaled);
end

% invert roundness
roundness_inv = 1./roundness;
% plot roundness...
roundness_sgf = sgolayfilt(roundness_inv(2:final_idx), sgf_deg, sgf_wind);
p7 = plot(tspan_scaled(2:final_idx), roundness_sgf, 'LineWidth', 1);
roundness_color = graphics_color('sky blue');
p7.Color = roundness_color;

% ... with error shade
std_wind = 21;
error_color = roundness_color.*0.7;
roundness_std = movstd(roundness_inv(2:final_idx),std_wind);
errorshade(tspan_scaled(2:final_idx), roundness_sgf - roundness_std, roundness_sgf + roundness_std, error_color)

if show_full_time
    axis([0,tspan_scaled(end), 0, 1]) 
else % show until final tstamp
    axis([0,tspan_scaled(tstamps(end)), 0, 1])
end

if plotleft
    set(gca,'ytick',[]) % some whitespace left?
else
    ylabel('Roundness',  'Interpreter','latex')
end

lgd = legend('$V_p$', '$V_q$', ...
    '$V_n$', 'Interpreter','latex');

yyaxis left
axis([0, 40.5, 0, 0.5])

%lgd.NumColumns = 2;
%set(gca,'xtick',[10 20 30 40])
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340/sqrt(2) 240/sqrt(2)]);
set(gca, 'fontname', 'Roman', 'FontSize', 9.0)
% uncomment to save
% exportgraphics(gca,'regionaldynamics_test10.pdf')

%% Boundary perturbation
% Update this to use function GrowthFactors.m
%error('Boundary perturbation section under construction!')

% estimate Lambda(k) from numerical boundary perturbation
J = 10;
tint = 5; % time intervals
sgf_wind = 21;   % savitzky-golay filter window
sgf_deg = 7;
lam_est = zeros(floor((numel(tspan)-2)/tint),J);
amp_k = zeros(numel(tspan)-1,J);
amp_k_sgf = zeros(numel(tspan)-1,J);

% find lambda estimates for each mode k = j
for j = 1:J    % j=0 should correspond to volume growth?
for i = 2:numel(tspan)
    amp_k(i-1,j) = BFTsave{i}(2,j+1);   % find amplitude evolution ...
    %Namp_k(i-1,j) = BFTsaveN{i}(2,j+1);     % only for testing, N, Q perts
    %Qamp_k(i-1,j) = BFTsaveQ{i}(2,j+1); 
end
amp_k_sgf(:,j) = sgolayfilt(amp_k(:,j), sgf_deg, sgf_wind); % ... and filter it

% if negative amplutides in filtering
if min(amp_k_sgf(:,j)) < 0
    warning('k = ' + string(j) + ' not filtered due to it resulting in negative amplitudes')
    amp_k_sgf(:,j) = amp_k(:,j);
end

% estimate lambda using the amplitude of mode j
lam_est(:,j) = log(amp_k_sgf(1+tint:tint:end, j) ./ amp_k_sgf(1:tint:end-tint, j)) ... 
                ./diff(tspan(2:tint:end)');
end

% moving standard deviation window
std_wind = 7;

i = 2;   % subfigure index
check_modes = [1 2 3 4];
% estimate 'analytical' Lambda(k) directly from numerical regional volumes

% get estimated regional perturbations
rqc = -1./j*rp.*(rp.^-j - rp.^j)./(rn.^-j - rn.^j).*rq./ ... 
    (rq.^2 - rn.^2).*(rn.^j./rq.^j - rq.^j./rn.^j);
rnc = rp./rn.*(rp.^-j - rp.^j)./(rn.^-j - rn.^j);
rqc(isnan(rqc)) = 0;    % remove NaN
rnc(isnan(rnc)) = 0;
oxypert = -rp.^(-j-1).*rq.^(j+1)*mu_prol.*rqc - ...  % oxygen perturbation effect
       rp.^(-j-1).*rn.^(j+1)*mu_death.*rnc;

LAM = mu_prol*0.5*(mu_death/mu_prol*rn.^2.*(1+j) + rq.^2.*(1+j) ...
                 + rp.^2.*(1-j))./rp.^2 - sigma*j.*(j.^2-1)./rp.^3 ...
                 + oxypert;

LAM(1) = []; % 'pop' the first reduntant entry (with rp = 0)
% get associated moving standard deviation for errorshade
LAMj_std = movstd(LAM,std_wind);
lam_estj_std = movstd(lam_est(:,j),std_wind);

% apply Savitzky-Golay filter on analytical lambda estimate
LAMj_sgf = sgolayfilt(LAM, sgf_deg, sgf_wind);
%LAMj_sgf = LAMj_sgf(half_sgf_wind:end);
lam_estj_sgf = sgolayfilt(lam_est(:,j), sgf_deg, sgf_wind);
%lam_estj_sgf = lam_estj_sgf(half_sgf_wind:end);

% plot for k = all j
figure(i)
colororder({'[0 0.4470 0.7410]', '[0.8500 0.3250 0.0980]', '[0 0.4470 0.7410]'})
plot(tspan(2:end-1), LAMj_sgf, 'lineWidth', 2)
hold on
error_color = [0 0.4470 0.7410];
errorshade(tspan(2:end-1), LAMj_sgf - LAMj_std, LAMj_sgf + LAMj_std, error_color)

% note: lam_est is evaluated in a simple manner at the left node
plot(tspan(2:tint:end-1), lam_estj_sgf, '--', 'lineWidth', 2)
error_color = [0.8500 0.3250 0.0980];
errorshade(tspan(2:tint:end-1), lam_estj_sgf' - lam_estj_std', lam_estj_sgf' + lam_estj_std', error_color)

title('k = ' + string(j), 'Interpreter','latex')
xlabel('Time [$\mu_{prol.}^{-1}$]', 'Interpreter','latex')
ylabel('$\Lambda(k;t)$ [$\rho_{prol.}$]', 'Interpreter','latex')
axis([0, tspan(end), min(lam_est(:,j))*1.0, max(lam_est(:,j))*1.2])

% yyaxis right
% p7 = plot(tspan(2:end), roundness(2:end), 'LineWidth', 1);
% p7.Color = graphics_color('sky blue');
% ylabel('Perimiter$^2/(4\pi*$Area$)$', 'Interpreter','latex')
% 
% legend('expected ', 'numerical', 'boundary roundness', 'Location','northwest', 'Interpreter','latex')

%sgtitle('Boundary perturbation dynamics, IC: k = ' + string(k), 'Interpreter','latex')
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 340/sqrt(2) 240/sqrt(2)]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)

i = i+1;    % next figure

