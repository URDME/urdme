function mexexec = mexmake_dlcm(propensity_file,varargin)
%MEXMAKE_DLCM Makefile for MEXDLCM.
%   MEXEXEC = MEXMAKE_DLCM(P) Makes the DLCM-solver with propensity
%   source file P, given as a relative path, and returns the
%   executable file as a function handle.
%
%   MEXEXEC = MEXMAKE_DLCM(P,...) accepts additional arguments as
%   property/value-pairs.
%
%   Property     Value/{Default}     Description
%   -----------------------------------------------------------------
%   mexname      string {'mexdlcm'}  Name of output mex-function, cf.
%                                    "mex -output <mexname>"
%
%   mexhash      scalar uint32 {0}   Model hashkey, cf.
%                                    "-DMEXHASH=<mexhash>"
%
%   define       string              Defines on the mex-format
%                                    '-D<define>', see MEX.
%   source       string              Source filename(s).
%   include      string              Include path.
%   link         string              Linker path.
%
%   See also DLCM, UDS, MEX.

% E. Blom 2024-11-25 (directly using mexmake_uds.m)

% reuse make from uds-solver...
mexmake_uds(propensity_file,varargin{:});

% ... but use proper mexexec
mexexec = str2func('mexdlcm');
