function URDMEtest
%URDMEtest Run all available URDME tests.

% S. Engblom 2019-11-14 (Minor revision)
% S. Engblom 2018-06-28

ok = 1;

% check that we are in the test-folder
path = mfilename('fullpath');
path = path(1:end-numel(mfilename)-1);
folder = pwd;
if ~strcmp(folder,path)
  warning('Tests should run from the test-folder.');
end

fprintf('\n\n*** Running all URDME tests... ***\n\n');
ftest = {'spin' 'spin_aem' 'spin_ssa' 'spin_uds' ...
         'spin_subdiffusion' ...
         'examples' ...
         'scrub' 'scrub_subdiffusion' 'scrub_neuron' 'scrub_DLCM'};
oktest = zeros(1,numel(ftest));
for i = 1:numel(ftest)
  oktest(i) = feval(ftest{i});
  fprintf('\n');
end

ok = all(oktest);
if ok
  fprintf('\n\n*** All URDME tests passed. ***\n\n');
else
  for i = find(~oktest)
    fprintf(sprintf('\n\n*** URDME test ''%s'' failed. ***\n',ftest{i}));
  end
  fprintf('\n');
end
