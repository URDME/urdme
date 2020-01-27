function [REAC,SPEC,RATE] = AMPA_NMDA_SYNAPSE(umod)
%   AMPA_NMDA_SYNAPSE builder.
%
%   Syntax [REAC,SPEC,RATE] = AMPA_NMDA_SYNAPSE or UMOD =
%   AMPA_NMDA_SYNAPSE(UMOD).

% S. Engblom 2019-12-06

% Kinetic scheme from:
% http://link.springer.com/article/10.1007/BF00961734

[REAC_AMPA,SPEC_AMPA,RATE_AMPA] = AMPA_SYNAPSE;
[REAC_NMDA,SPEC_NMDA,RATE_NMDA] = NMDA_SYNAPSE;
REAC = [REAC_AMPA REAC_NMDA];
SPEC = [SPEC_AMPA SPEC_NMDA];
% shared transmitter "T" passed in ldata, here removed from AMPA
% to avoid the duplicate in NMDA:
RATE = [RATE_AMPA(1:end-2) RATE_NMDA];

if nargin > 0
  if nargout <= 1
    umod = rparse(umod,REAC,SPEC,RATE,'synapse/lib/AMPA_NMDAprop.c');
    % additional small library included
    rparse_include('synapse/lib/AMPA_NMDAprop.c','#include "NMDAfun.h"');
    umod.makeargs = {'source' 'synapse/lib/NMDAfun.c'};
    REAC = umod;
  else
    error('Unknown syntax.');
  end
end
