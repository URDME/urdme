function R = mexrhs
%MEXRHS Right-hand side evaluation of URDME propensities.
%   R = MEXRHS(MEXHASH,tspan,X,size(N,2),VOL,LDATA,GDATA, ...
%     LDATA_TIME,GDATA_TIME,SD,K,I,S)
%   returns an (Mreactions*Ncells)-by-NUMEL(tspan) matrix R such that
%   R(:,k) contains the values of all Mreactions propensities in all
%   Ncells voxels at time tspan(k). Let R_ =
%   reshape(R(:,k),Mreactions,Ncells). Then R_(i,j) contains the value
%   of propensity #i in voxel j.
%
%   The input state X is an Mspecies-by-Ncells matrix and the format
%   for the other inputs are described in URDME and URDME_VALIDATE.
%
%   MEXRHS also supports a syntax with multiple replicas for which
%   size(R,3) = size(X,3), the number of different replicas.
%
%   Note: MEXRHS is intended for internal use and no error-checking is
%   performed.
%
%   See also MEXJAC, UDS, URDME, URDME_VALIDATE.

% S. Engblom 2020-02-21

error('This file should not be called directly.');
