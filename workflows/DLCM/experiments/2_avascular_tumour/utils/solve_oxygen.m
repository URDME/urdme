function Oxy = solve_oxygen(iOLa, M, Oxy, U, cons, extdof, cutoff_bdof, cutoff_die, dtau, dt)
% SOLVE_OXYGEN Solves the oxygen problem given no necrotic consumption 
% in the tumor model
% Assumes fixed tumor domain, given by U. Finds the oxygen solution by
% time stepping the initial guess, Oxy,
% using Euler backward with explicit source term for the time-dependent
% heat equation.
%      
%   input:  iOLa    - LU-factorization structs for computational domain
%           M       - Mass matrix
%           Oxy     - Initial guess, oxygen field
%           U       - Cell density field
%           cons    - Oxygen consumption rate
%           extdof  - Domain boundary node indices
%           cutoff_bdof - Cell density cutoff in PDE solver
%           cutoff_die  - Oxygen threshold below which cells die
%           dtau    - Oxygen solver timestep (pseudo-time)
%           dt      - PDE solver timestep
%   output: Oxy     - Approximated oxygen solution at t+dt
%

% E. Blom 2023-01-01

% for radially symmetric tumors, 
% this returns values that agree well with the analytical oxygen relations
% kappa_die  ~
%    cons*(0.5*rp^2*log(rp)-0.5*rn^2*log(rn) + 0.25*(rn^2 - rp^2)) + 1
% kappa_prol ~ 
%   cons*(0.5*rp^2*log(rp)-0.5*rn^2*log(rq) + 0.25*(rq^2 - rp^2)) + 1

%dtau= 1e-3;     %delta tau/epsilon 
%MinvA = M\A;               % explicit scheme
%MinvA = (M + dtau*A)\M;     % Euler backward scheme

Oxy_current = Oxy;
Oxy_new = zeros(numel(Oxy),1);

% find solution iteratively using Backwards Euler scheme
for t = 1:round(dt/dtau) + 1
    
    % consider the source term explicitly
    eta_star = FindEtaStar(Oxy_current, U, cons, extdof, cutoff_bdof, cutoff_die);
    
    %
    % Backward Euler
    Oxy_new(iOLa.q) = iOLa.U\(iOLa.L\(iOLa.R(:,iOLa.p)\M*(Oxy_current - dtau.*eta_star)));
    Oxy_diff = Oxy_new - Oxy_current;
    err(t) = sqrt(Oxy_diff'*M*Oxy_diff);         % local error
    %differ(t) = sum(abs(Oxy_new - Oxy)); % difference between intitial guess
    Oxy_current = Oxy_new;
    
    % check stationary solution residual
    % must adapt eta_star at the boundary for this
    %eta_star_expl = eta_star;
    %eta_star_expl(extdof) = 1; % rest is the same
    
    %stationary_residual = A*Oxy_current - M*eta_star_expl;
    %stationary_residual(extdof) = 0;    % exclude boundary as eq. breaks there
    %residual_norm(t) = sum(abs(stationary_residual));
    
end

Oxy = Oxy_new; % return the solution
% indication of bad convergence
if max(err)/dtau > 1000
    warning('large local errors for the iterative oxygen solver')
end

% finds the explicit source term that depends on Oxy_current
function eta_star = FindEtaStar(oxy, u, cons, extdof, cutoff_bdof, cutoff_die)

    % evaluate necrotic dofs
    % Why do not all voxels below cutoff become ndofs here...
    ndof = find(u & (oxy < cutoff_die)); % 'necrotic' dof
    
    eta_star = cons.*u./u; % consumption sink term
    
    eta_star(ndof) = 0.*eta_star(ndof); % necrotic cells consume no oxygen
    % constant cell density during iterative solver
    eta_star(u <= cutoff_bdof) = 0;    % no consumption outside artificial tumor boundary
    eta_star(extdof) = -1;% -1 for BE and 1 for explicit methods.

end
end