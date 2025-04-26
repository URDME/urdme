function sol = RegionalDynamics2D(par,Vp0,tspan)
%REGIONALDYNAMICS2D Solution of 2D radially symmetric tumor.
%   Input:
%     par.mu_prol     - rate of proliferation
%     par.mu_death    - rate of cell death
%     par.kappa_prol  - oxygen threshold for proliferation
%     par.kappa_death - oxygen threshold for cell death
%     par.lambda      - rate of oxygen consumption per cell unit density
%     par.c_out       - oxygen boundary condition at VR
%     Vp0             - initial tumor volume
%     tspan           - simulation times output
%
%   Output:
%     sol.Vp    - total volume
%     sol.Vq    - volume quiescsent + necrotic region
%     sol.Vn    - volume necrotic region
%     sol.discr - stability discriminant
%     sol.tspan - time steps used in the solver

% S. Engblom 2023-05-12 (Revision)
% E. Blom 2023-01-01

assert(par.mu_prol == 1,'Incorrect non-dimensionalization (mu_prol)');
assert(par.c_out == 1,'Incorrect non-dimensionalization (c_out)');

% non-dimensional parameters
K_prol = 4*(1-par.kappa_prol)/par.lambda;
K_death = 4*(1-par.kappa_death)/par.lambda;

% solve as an ODE
[tspan,rp2] = ode15s(@l_RHS,tspan,Vp0/pi, ...
                     odeset('RelTol',1e-5,'AbsTol',1e-8), ...
                     K_prol,K_death,par.mu_death);

% pick out radii and postprocess
[rn2,rq2] = l_rads2(rp2,K_prol,K_death);
sol.Vp = pi*rp2;
sol.Vq = pi*rq2;
sol.Vn = pi*rn2;
sol.discr = l_discriminant(rn2,rq2,rp2,par.mu_death);
sol.tspan = tspan;

end

% ----------------------------------------------------------------------
function F = l_RHS(t,rp2,K_prol,K_death,mu_death)
%L_RHS ODE RHS equation.

[rn2,rq2] = l_rads2(rp2,K_prol,K_death);
F = -(mu_death*rn2+rq2-rp2);

end
% ----------------------------------------------------------------------
function D = l_discriminant(rn2,rq2,rp2,mu_death)
%L_DISCRIMINANT Discriminant.
%   Note: assumes all r-vectors to have the same sizes.

D = 1-log(rp2)./log(rn2).*(mu_death+log(rq2./rn2)./(rq2-rn2).*rq2);
  
end
% ----------------------------------------------------------------------
function [rn2,rq2] = l_rads2(rp2,K_prol,K_death)
%L_RADS2 Determine rn2 and rq2 from rp2 and parameters K_prol, K_death.

D = @(r2)r2.*log(r2)-r2;
D0 = D(rp2);
D2 = -K_death-D0;
D1 = -K_prol-D0;

% rn2 = rn^2 first, avoid complex values
ix = 0 <= D2 & D2 <= 1;
rn2 = zeros(size(rp2)); % ensure 0 <= rn2

% formula:
arg = -D2(ix)*exp(-1);
rn2(ix) = exp(lambertw(-1,arg)+1);
rn2 = real(rn2); % clear roundoff errors

% Identity: exp(W(z)) = z/W(z)
%
% rq2 = -rn2*W(-1/rn2*exp(-D1/rn2))
%
% exp(-D1/rn2) = exp(-D1*exp(-1)*W(X)/X) =
% exp(D1/D2*W(X)) = rn2^(D1/D2)*exp(-D1/D2)
%
% -1/rn2*exp(-D1/rn2) = -rn2^(D1/D2-1)*exp(-D1/D2) = 
% -exp(-1)*Y^(D1/D2-1), Y := rn2*exp(-1)
% -exp(1)*Y*W(-exp(-1)*Y^(D1/D2-1))
%
% = -exp(1)*Y*[-exp(-1)*Y^(D1/D2-1)] * 
%   W(-exp(-1)*Y^(D1/D2-1))/[-exp(-1)*Y^(D1/D2-1)]
% = Y^D1/D2 * exp(-W(-exp(-1)*Y^(D1/D2-1)))
% = exp(1)*Y * exp(-1)*Y^(D1/D2-1) * exp(-W(-exp(-1)*Y^(D1/D2-1)))
% = Y*Z*exp(1-W(-Z)), Z := Y^(D1/D2-1)*exp(-1)

% rq = rq^2, again avoid complex values
ix = rn2 > 0;
rq2 = zeros(size(rp2));
rq2(~ix) = max(D1(~ix),0); % ensure 0 <= rq2

% use expansion for small values of rn2:
ixa = ix & rn2 < D1/100;
ixb = ix & rn2 >= D1/100;

% use the series expansion lambertw(-1,-x) = log(x)-log(-log(x))+R,
% where R ~ log(-log(x))/log(x) ~ <small>...
arg = rn2(ixa).*log(rn2(ixa))+D1(ixa);
rq2(ixa) = D1(ixa)+rn2(ixa).*log(arg);

% second formula:
%rq2(ixb) = -rn2(ixb).*lambertw(-1,-1./rn2(ixb).*exp(-D1(ixb)./rn2(ixb)));
% a bit prettier...
%arg = rn2(ixb)*exp(-1);
%rq2(ixb) = -exp(1)*arg.*lambertw(-1,-exp(-1)*arg.^(D1(ixb)./D2(ixb)-1));
% ...and even somewhat prettier still:
arg1 = rn2(ixb)*exp(-1);
arg2 = arg1.^(D1(ixb)./D2(ixb)-1)*exp(-1);
rq2(ixb) = arg1.*arg2.*exp(1-lambertw(-1,-arg2));
rq2 = real(rq2); % clear roundoff errors

% finalize
rq2 = min(rp2,rq2); % ensure rq <= rp

end
% ----------------------------------------------------------------------
  
% $$$ %#### STATUS ####
% $$$ % The Euler time-Stepping becomes highly oscillatory after Vp > 1
% $$$ % might be a problem for some parameter regimes also within feasible
% $$$ % domain?
% $$$ % 
% $$$ % NOTE: the scaling VR is important when comparing solutions
% $$$ % using radial scaling R = 1.
% $$$ 
% $$$ % Add tester to the solver that checks
% $$$ % whether or not the solution values satisfy the
% $$$ % governing equations and so on
% $$$ %####
% $$$ 
% $$$ % The idea is that the solution is unique within the physically
% $$$ % feasible domain, i.e, Vn <= Vq <= Vp.
% $$$ 
% $$$ 
% $$$ VR = pi;    % keeps scaling assumption that R = 1
% $$$ 
% $$$ % define constants of nonlinear regional relations 
% $$$ c1 = 4*pi*(par.kappa_prol - par.kappa_die)/par.lambda;
% $$$ c2 = 4*pi*(par.kappa_die - par.c_out)/par.lambda;
% $$$ c3 = 4*pi*(par.kappa_prol - par.c_out)/par.lambda;
% $$$ 
% $$$ % define nonlinear regional relations
% $$$ v1 = @(vq, vn)abs(vq + vn.*(log(vn./vq)-1) - c1);                % relates Vq and Vn
% $$$ v2 = @(vp, vn)abs(vp.*(log(vp/VR) - 1) - vn.*(log(vn/VR) - 1) - c2);   % relates Vp and Vn
% $$$ v3 = @(vp, vq)abs(vp.*(log(vp/VR) - 1) + vq - c3);   % relates Vp and Vq for Vn == 0
% $$$ 
% $$$ % define discriminants
% $$$ discr1 = @(vp, vq, vn)(par.r_prol - log(vp/VR)./log(vn/VR)...
% $$$     .*(par.r_die + par.r_prol.*vq.*log(vq./vn)./(vq - vn)));  % volume discriminant
% $$$ discr2 = @(vp, vq)(par.r_prol*(1 + log(vp/VR)));
% $$$ 
% $$$ % volume increase
% $$$ vvel = @(vp, vq, vn)(-par.r_die.*vn - par.r_prol.*vq + par.r_prol.*vp);
% $$$ 
% $$$ %% Plot dynamics for an interval [0,T]
% $$$ % test simple Euler forward for simulating the region dynamics
% $$$ T = Tend;
% $$$ N = 3000;
% $$$ dt = T/N;
% $$$ tspan = 0:dt:T;
% $$$ 
% $$$ discriminant = 0; % reset variable 'discriminant'
% $$$ % initial condition
% $$$ Vp(1) =  Vp0;
% $$$ for i = 1:N
% $$$     
% $$$     % Solve AE for Vq, Vn given Vp
% $$$     [Vq(i), Vn(i)] = FindQuiescentAndNecroticVolumes(Vp(i), v2, v1, v3);
% $$$     
% $$$     Vp(i+1) = Vp(i) + vvel(Vp(i), Vq(i), Vn(i))*dt;    % Euler forward
% $$$     
% $$$     % evaluate radial discriminant
% $$$     if Vn(i) > 0
% $$$         discriminant(i) = discr1(Vp(i), Vq(i), Vn(i));
% $$$     elseif Vq(i) > 0
% $$$         discriminant(i) = discr2(Vp(i), Vq(i));
% $$$     else
% $$$         discriminant(i) = par.r_prol; % just a positive constant...
% $$$     end
% $$$     
% $$$     if Vp(i+1) > VR
% $$$         warning('Simulation aborted. Tumor reached oxygen source.');
% $$$         tspan = tspan(1:i+1);   % truncate tspan vector
% $$$         break;
% $$$     end 
% $$$ end
% $$$ 
% $$$ % truncate sizes
% $$$ Vp = Vp(1:end-1);
% $$$ tspan = tspan(1:end-1);
% $$$ %i = i+1;
% $$$ %Vn(i) = fzero(@(vn)(v2(Vp(i),vn)), [0.001 Vp(i)]);
% $$$ %Vq(i) = fzero(@(vq)(v1(vq,Vn(i))), [Vn(i) Vp(i)]);
% $$$ 
% $$$ sol = struct('Vp', {Vp}, 'Vq', {Vq}, 'Vn', {Vn}, ...
% $$$     'radial_discriminant', {discriminant}, 'tspan', {tspan});
% $$$ 
% $$$ return;
% $$$ 
% $$$ %% Plots used for debugging
% $$$ % time intervals for plotting markers
% $$$ tint = round(N/10);
% $$$ t_idx = sort([2:tint:numel(tspan)]);
% $$$ 
% $$$ % plot the regional characteristics
% $$$ figure(1)
% $$$ p1 = plot(tspan(t_idx), Vp(t_idx), '*', 'LineWidth', 1);
% $$$ hold on
% $$$ p2 = plot(tspan(t_idx), Vq(t_idx), '.', 'Marker', '.', 'MarkerSize', 15);
% $$$ plot(tspan(t_idx), Vn(t_idx), 'kx', 'LineWidth', 1)
% $$$ p1.Color = graphics_color('vermillion');
% $$$ p2.Color = graphics_color('bluish green');
% $$$ p3 = plot(tspan, Vp, '-', 'LineWidth', 1);
% $$$ p4 = plot(tspan, Vq, '-', 'LineWidth', 1);
% $$$ plot(tspan, Vn, 'k-', 'LineWidth', 1);
% $$$ p3.Color = graphics_color('vermillion');
% $$$ p4.Color = graphics_color('bluish green');
% $$$ %plot(tspan(1:end-1), discriminant, '--', 'LineWidth', 2)
% $$$ %plot(tspan(1:end-1), sonr + 0.*tspan(1:end-1), 'k--')
% $$$ xlabel('Time [$\rho_{prol.}^{-1}$]', 'Interpreter','latex')
% $$$ ylabel('Volume [$V_R$]', 'Interpreter','latex')
% $$$ %axis( [0,tspan(end), 0,0.4])
% $$$ 
% $$$ yyaxis right
% $$$ plot(tspan(t_idx), 0.5*vvel(Vp(t_idx), Vq(t_idx), Vn(t_idx))./Vp(t_idx), '--', 'LineWidth', 2)
% $$$ ylabel('$\Delta_{\theta} [\mu_{prol.}]$', 'Interpreter', 'Latex')
% $$$ 
% $$$ legend('$V_p$', '$V_q$', '$V_n$', 'Interpreter', 'Latex')
% $$$ 
% $$$ 
% $$$ set(gcf,'PaperPositionMode','auto');
% $$$ set(gcf,'Position',[100 100 340/(sqrt(2)) 240/(sqrt(2))]);
% $$$ set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
% $$$ % uncomment to save
% $$$ % exportgraphics(gca,'regionaldynamics_ex.pdf')
% $$$ 
% $$$ 
% $$$ % plot analytical perturbation amplification
% $$$ modes = [1 3 9];
% $$$ warning('sigma manually input!')
% $$$ Lambda = PerturbationAmplificationFactors(sqrt(Vp/pi), sqrt(Vq/pi), ...
% $$$                                 sqrt(Vn/pi), 0.00000, par.r_prol, par.r_die, modes);
% $$$ figure(2)
% $$$ % plot markers
% $$$ plot(tspan(t_idx), Lambda(modes(1),t_idx), 'o', 'LineWidth', 1)
% $$$ hold on
% $$$ plot(tspan(t_idx), Lambda(modes(2),t_idx), '^', 'LineWidth', 1)
% $$$ plot(tspan(t_idx), Lambda(modes(3),t_idx), 'square', 'LineWidth', 1)
% $$$ %plot(tspan(2:tint:end-1), Lambda(modes(4),2:tint:end), '-pentagram', 'LineWidth', 1, 'MarkerSize', 10)
% $$$ %hold on; 
% $$$ % plot positive term for k=1 (markers)
% $$$  plot(tspan(t_idx), par.r_prol*0.5*(par.r_die/par.r_prol*Vn(t_idx).*(1+1) + Vq(t_idx).*(1+1) ...
% $$$                           + Vp(t_idx).*(1-1))./Vp(t_idx), 'k+', 'MarkerSize', 10)
% $$$  % % plot negative terms for k=1 
% $$$  % yyaxis right
% $$$  plot(tspan(t_idx), -(Lambda(1,t_idx) - par.r_prol*0.5*(par.r_die/par.r_prol*Vn(t_idx).*(1+1) + Vq(t_idx).*(1+1) ...
% $$$                           + Vp(t_idx).*(1-1))./Vp(t_idx)), 'k_', 'MarkerSize', 10)
% $$$  % plot full lines
% $$$  plot(tspan, Lambda(modes(1),:), '-', 'LineWidth', 1,  'color', [0 0.4470 0.7410])
% $$$  hold on
% $$$  plot(tspan, Lambda(modes(2),:), '-', 'LineWidth', 1,  'color', [0.8500 0.3250 0.0980])
% $$$  plot(tspan, Lambda(modes(3),:), '-', 'LineWidth', 1,  'color', [0.9290 0.6940 0.1250])
% $$$  plot(tspan, par.r_prol*0.5*(par.r_die/par.r_prol*Vn.*(1+1) + Vq.*(1+1) ...
% $$$                           + Vp.*(1-1))./Vp, 'k:')
% $$$  plot(tspan, -(Lambda(1,:) - r_prol*0.5*(par.r_die/par.r_prol*Vn.*(1+1) + Vq.*(1+1) ...
% $$$                           + Vp.*(1-1))./Vp), 'k:')
% $$$ 
% $$$  xlabel('Time [$\rho_{prol.}^{-1}$]', 'Interpreter','latex')
% $$$  ylabel('$\Lambda(k)$', 'Interpreter','latex')
% $$$ % 
% $$$  axis([0,Tend, 0, 1.1])
% $$$  hold on
% $$$  
% $$$ % axis([0,tspan(end), -1 0])
% $$$ % ylabel('negative contribution, $k=1$', 'Interpreter','latex')
% $$$ legend('$k=1$', '$k=3$', '$k=9$', 'Interpreter','latex')
% $$$ 
% $$$ set(gcf,'PaperPositionMode','auto');
% $$$ set(gcf,'Position',[100 100 340/(sqrt(2)) 240/(sqrt(2))]);
% $$$ set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
% $$$ % uncomment to save
% $$$ % exportgraphics(gca,'lambdadynamics_ex.pdf')
% $$$ 
% $$$ end
