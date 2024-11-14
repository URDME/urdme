function JAC = mexjac
%MEXJAC Jacobian for propensities in URDME.
%   JAC = MEXJAC(MEXHASH,T,X,size(N,2),G,VOL,LDATA,GDATA, ...
%     LDATA_TIME,GDATA_TIME,SD,K,I,S)
%   returns a reaction Jacobian JAC, an
%   (Mreactions*Ncells)-by-(Mspecies*Ncells) sparse matrix which is
%   the Jacobian of MEXRHS with the same arguments (except that MEXJAC
%   requires also the sparse dependency matrix G as input).
%
%   Notes: MEXJAC only handles a scalar time T. Inline propensities
%   are always supported while for compiled propensities, the
%   propensity file has to be augmented with a GET_jacobian()-function
%   and comiled with JAC_ #define'd. MEXJAC is intended for internal
%   use and no error-checking is performed.
%
%   See also MEXRHS, UDS.

% S. Engblom 2024-06-10
% S. Engblom 2020-02-21

error('This file should not be called directly.');
