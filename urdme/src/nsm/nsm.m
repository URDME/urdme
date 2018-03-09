function nsm
%NSM URDME Next Subvolume Method solver.
%   The URDME NSM solver is the default simulation engine in URDME and
%   has been present since URDME 1.0. It implements the "Next
%   Subvolume Method" [1], but features a detailed dependency graph G
%   including both diffusion and reaction events. Hence the URDME NSM
%   has the ability to use computed rates in an optimal way.
%
%   The NSM solver supports three report levels: 0 (no report), 1
%   (progress report), 2 (progress report + event count).
%
%   Property    Value/{Default}           Description
%   -----------------------------------------------------------------------
%   K, I, S     Matrices                  Inline propensities.
%
%   The NSM solver supports the 'inline' propensity syntax, where
%   UMOD.solverargs = {'K' K 'I' I 'S' S}. This is an efficient way of
%   defining elementary reactions and is parsed as follows:
%
%   -K, the inline propensity coefficient matrix. K(:,i) gives the
%   coefficients to reaction #i. Where i is in [1,M1]. The size of K
%   should be 3-by-M1.
%
%   -I, the inline propensity index matrix. I(:,i) gives the indices
%   used by reaction #i where i is in [1,M1]. The size of I should be
%   3-by-M1.
%
%   -S (optional), the inline propensity subdomain matrix. S(:,i)
%   lists all subdomains in which the corresponding inline propensity
%   is off. S should be a sparse matrix of size <anything>-by-M1.
%
%   The inline propensities are defined by three vectors. [k1,k2,k3],
%   [i,j,k] and [s1,s2...sn] taken from the columns of K, I and S,
%   respectively.
%
%   If i == j the result is
%     k1*x[i]*(x[i]-1)/(2*vol) + k2*x[k] + k3*vol,
%   else if i ~= j,
%     k1*x[i]*x[j]/vol + k2*x[k] + k3*vol.
%   Note the scaling with 2 for the case i == j.
%
%   An inline propensity is automatically off (the result is zero) in
%   subdomain zero. Additionally it is off in all subdomains listed in
%   [s1,s2...sn].
%
%   Example:
%     % the reversible transformation X+Y <--> X+Z 
%     % under 1D periodic transport of X
%     clear umod;
%
%     % periodic transport operator over 5 voxels
%     D = [-1  1  0  0  0; ...
%           0 -1  1  0  0; ...
%           0  0 -1  1  0; ...
%           0  0  0 -1  1; ...
%           1  0  0  0 -1];
%     D3 = sparse(15,15); % no transport for Y and Z
%     D3(1:3:15,1:3:15) = D;
%     umod.D = D3';
%
%     % reaction topology
%     umod.N = sparse([ 0  0; ...
%                      -1  1; ...
%                       1 -1]);
%     umod.vol = ones(1,5);
%     umod.sd = ones(1,5);
%
%     % inline propensities
%     K = [10 0  0; ... % bimolecular
%          0  2  0]';   % decay
%     I = [1 2 1; ...
%          1 1 3]';
%     umod.solverargs = {'K' K 'I' I};
%
%     % dependency graph
%     umod.G = sparse([1 1 0 1 1; ...
%                      0 0 1 1 1]);
%
%     % initial state u0 and time interval
%     umod.u0 = [1 0 0 0 0; ...
%                1 1 1 1 1; ...
%                0 0 0 0 0];
%     umod.tspan = 0:0.5:25;
%
%     % compile and solve...
%     umod = urdme(umod,'solver','nsm');
%
%     % plot Z in 1st voxel only
%     figure, stairs(umod.tspan,umod.U(3,:));
%     ylim([-0.1 1.1]);
%     % (the X-molecule is expected to complete one round in about 5
%     % units of time, giving rise to a quasi-periodic behavior)
%
%   See also UDS, RPARSE_INLINE.
%
%   References:
%     [1] D. Fange and J. Elf: "Noise-Induced Min Phenotypes in
%     E. coli", PLoS Comput. Biol. 2(6):637--648 (2006).
                
% S. Engblom 2017-02-22

error('This file should not be called directly.');
