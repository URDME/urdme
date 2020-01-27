function [REAC,SPEC,RATE] = HHKNa_CHANNEL(umod)
%   HHKNa_CHANNEL builder.
%
%   Syntax [REAC,SPEC,RATE] = HHKNa_CHANNEL or UMOD =
%   HHKNa_CHANNEL(UMOD).

% S. Engblom 2019-12-06

[REACK,SPECK,RATEK] = HHK_CHANNEL;
[REACNa,SPECNa,RATENa] = HHNa_CHANNEL;
REAC = [REACK REACNa];
SPEC = [SPECK SPECNa];
% shared voltage "v" passed in ldata, here removed from K to avoid the
% duplicate in Na:
RATE = [RATEK(1:end-2) RATENa];

if nargin > 0
  if nargout <= 1
    umod = rparse(umod,REAC,SPEC,RATE,'channel/lib/HHKNa_prop.c');

    % additional small libraries included
    rparse_include('channel/lib/HHKNa_prop.c','#include "HHNafun.h"');
    rparse_include('channel/lib/HHKNa_prop.c','#include "HHKfun.h"');
    umod.makeargs = {'source' ...
                     {{'channel/lib/HHKfun.c' 'channel/lib/HHNafun.c'}}};
    REAC = umod;
  else
    error('Unknown syntax.');
  end
end
