function [REAC,SPEC,RATE] = HHK_CHANNEL
%   HHK_CHANNEL builder.
%
%   Syntax [REAC,SPEC,RATE] = HHK_CHANNEL or UMOD = HHK_CHANNEL(UMOD).

% S. Engblom 2019-12-02 (Revision, previously ionprops.m)
% S. Engblom 2018-06-27 (Revision)
% A. Senek 2017-05-31

RATE = {'v' 'ldata'}; % (note: input via ldata)

SPEC = {'K1' 'K2' 'K3' 'K4' 'K5'};

% Supporting functions (HHKfun.h and HHKfun.c)
REAC = {'K1 > 4 * HHK_alpha1(v) * K1 > K2' ...
        'K2 > 1 * HHK_beta1(v)  * K2 > K1' ...
        'K2 > 3 * HHK_alpha1(v) * K2 > K3' ...
        'K3 > 2 * HHK_beta1(v)  * K3 > K2' ...
        'K3 > 2 * HHK_alpha1(v) * K3 > K4' ...
        'K4 > 3 * HHK_beta1(v)  * K4 > K3' ...
        'K4 > 1 * HHK_alpha1(v) * K4 > K5' ...
        'K5 > 5 * HHK_beta1(v)  * K5 > K4' ...
       };

if nargin > 0
  if nargout <= 1
    umod = rparse(REAC,SPEC,RATE,'channel/lib/HHK_prop.c');
    % additional small library included
    rparse_include('channel/lib/HHKNa_prop.c','#include "HHKfun.h"');
    umod.makeargs = {'source' 'channel/lib/HHKfun.c'};
    REAC = umod;
  else
    error('Unknown syntax.');
  end
end
