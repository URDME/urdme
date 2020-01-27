function [gamma,A,p] = subdata(phi_)
%SUBDATA Load pre-simulated subdiffusion data.
%   [GAMMA,A,P] = SUBDATA(PHI), for fraction of occupied space PHI,
%   loads pre-simulated coarse-grained subdiffusion data and computes
%   the internal states scalar diffusion coefficients GAMMA, the
%   internal states transfer matrix A, and a steady-state probability
%   distribution P over the internal states.

% S. Engblom 2017-05-05 (Revision)
% S. Engblom 2017-03-21 (with data from L. Meinecke)

% load pre-simulated data
load data/CoarseGrainedSD

% find phi
[~,i] = min(abs(phi-phi_));
if abs(phi(i)-phi_) > 0.01*abs(phi_)
  warning('Close match to value of phi seems to be missing.');
end

% single out column vectors
p = P(:,i);
gamma = Gamma(:,i);
Mu = Mu(:,i);

% eq. (3.3) \ref{eq:diffintst} (with gamma0 = 1)
theta = gamma;

% eq. (3.4) \ref{eq:mudef} (with f --> p)
mu = p.*theta;
mu = mu/sum(mu);

% eq. (3.5) \ref{eq:Adef}
A = mu*theta';
A = A-diag(theta);
% (check: sum(A,1) = 0)
