function uds
%UDS URDME Deterministic Solver.
%   The UDS solves the URDME model in the sense of a linearized moment
%   equation (see [1]) using Matlab's ODE-solvers. This is useful for
%   judging stochastic effects and evaluate deterministic
%   approximations.
%
%   Property    Value/{Default}           Description
%   -----------------------------------------------------------------------
%   odesolv     ODE-solver {@ode23s}      ODE-solver. See ODESET.
%   odeopts     ODE-solver options
%               {odeset('RelTol',1e-4,
%                'AbsTol',0.1)}
%
%   jacobian    {0} | 1                   Use exact Jacobian.
%
%   finish      {''}                      Play sound when done. Portable 
%                                         choices include 'handel', 'gong' 
%                                         and 'splat'. See AUDIOVIDEO.
%
%   Using an exact Jacobian requires either that all propensities are
%   inline or that the Jacobian is prepared in the propensity.c-file
%   and compiled and linked. See MEXMAKE_UDS.
%
%   Example:
%     % a basic bimolecular model in 10 volumes
%     clear umod;
%
%     % diffusion with periodic boundaries
%     n = 10;
%     h = 1/(n-1);
%     B = repmat([1 -2 1]/h^2,n,1);
%     D = spdiags(B,[-1 0 1],n,n);
%     D(end,1) = 1/h^2;
%     D(1,end) = 1/h^2;
%     umod.D = kron(D,speye(2));
%
%     % scaling
%     umod.vol = repmat(h,1,n);
%     umod.sd = ones(1,n);
%
%     % bimolecular model
%     umod = rparse_inline(umod, ...
%              {'@ > k > X','@ > k > Y','X+Y > mu > @'}, ...
%              {'X' 'Y'},{'k' 1e2 'mu' 1e-2});
%
%     % initial state u0 and time interval
%     umod.u0 = 50*ones(2,n);
%     umod.tspan = linspace(0,10,100);
%
%     % compile and solve once...
%     umod = urdme(umod,'solver','uds','solverargs',{{'jacobian' 1}});
%     U1 = umod.U;
%
%     % ...solve with NSM
%     umod.solverargs = {};
%     umod = urdme(umod,'solver','nsm');
%     U2 = umod.U;
%
%     % plot total number of molecules
%     figure, plot(umod.tspan,sum(U1),'b'), hold on,
%     stairs(umod.tspan,sum(U2),'r');
%
%   See also RPARSE_INLINE, SSA, NSM.
%
%   References:
%     [1] S. Engblom: "Computing the Moments of High Dimensional
%     Solutions of the Master Equation", 
%     Appl. Math. Comput. 18(2):498--515 (2006).

% S. Engblom 2020-02-21 (Revision, Jacobian)
% S. Engblom 2019-11-27 (Revision, inline propensities)
% S. Engblom 2017-02-28

error('This file should not be called directly.');
