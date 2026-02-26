% STARTUP startup for the DLCM workflow.

% S. Engblom 2026-02-10 (Minor revision for U15)
% E. Blom 2025-11-14

% check if few stenglib functions exist, otherwise issue a warning
funstr = {'tprod' 'tsum' 'fsetop' 'fsparse'};
% function use-cases:
% * tprod used in mesh2dual
% * tsum used in mesh2dual
% * fsetop used in mesh2dual and map inside mexdlcm
% * fsparse used in mexdlcm and for model initial conditions
for s = 1:numel(funstr)
  if exist(funstr{s}) ~= 3
    warning(sprintf(['Function %s not found on path. ' ...
                     'Add path to stenglib for full ' ...
                     'DLCM functionality.'],funstr{s}));
  end
end
clear funstr s
