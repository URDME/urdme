function JAC = mexjac
%MEXJAC Jacobian for propensities in URDME.
%   JAC = MEXJAC(T,X,size(N,2),G,VOL,LDATA,GDATA,SD,K,I,S) returns a
%   reaction Jacobian JAC, an (Mreactions*Ncells)-by-(Mspecies*Ncells)
%   sparse matrix which is the Jacobian of MEXRHS with the same
%   arguments (except that MEXJAC requires also the sparse dependency
%   matrix G).
%
%   Note: MEXJAC only handles a scalar time T and only inline
%   propensities are currently supported. MEXJAC is intended for
%   internal use and no error-checking is performed.
%
%   See also MEXRHS, UDS.

% S. Engblom 2020-02-21

error('This file should not be called directly.');
