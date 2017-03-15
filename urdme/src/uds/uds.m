function uds
%UDS URDME Deterministic Solver.
%   The UDS solves the URDME model in the sense of a linearized moment
%   equation (see [1]) using Matlab's ODE-solvers.
%
%   Property    Value/{Default}           Description
%   -----------------------------------------------------------------------
%   odesolv     ODE-solver {@ode23s}      ODE-solver. See ODESET.
%   odeopts     ODE-solver options
%               {odeset('RelTol',1e-4,
%                'AbsTol',0.1)}
%
%   K, I, S     Matrices                  Inline propensities, see NSM.
%
%   finish      {''}                      Play sound when done. Portable 
%                                         choices include 'handel', 'gong' 
%                                         and 'splat'. See AUDIOVIDEO. 
%
%   See also RPARSE_INLINE, SSA, NSM.
%
%   References:
%     [1] S. Engblom: "Computing the Moments of High Dimensional
%     Solutions of the Master Equation"
%     Appl. Math. Comput. 18(2):498--515 (2006).

% S. Engblom 2017-02-28

error('This file should not be called directly.');
