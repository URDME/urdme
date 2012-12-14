%MINSWEEP.  Example of a parameter sweep based on the MINCDE model.
%
%   Example:
%     minsweep          % wait until done!
%
%     minsweep_plot(1); % inspect the first
%     minsweep_allplots
%
%   See also COLI_MODEL, HUANG, MINSWEEP_ALLPLOTS.

% P. Bauer   2012-10-26 (Update to Comsol 4.x)
% S. Engblom 2011-07-14 (Revision)
% S. Engblom and A. Hellander 2011-06-01

% specify comsol version - '3.5' or '4.x'
comsol ='4.x';

% parameter space to sweep over
%Nval = 30;
Nval=4;

xsep = linspace(0,4.5e-6,Nval+1);

xsep(end) = []; % (avoid creating two distinct bacteria)

if(~exist('results/','file')>0),system('mkdir results');end
save results/info.mat xsep % convenient

for i = 1:Nval
  if(file_exists(sprintf('results/out%d.mat',i)))
    fprintf('output file %s found, skipping\n',sprintf('results/out%d.mat',i));
    continue;
  end
  % generate two merged E. colis with separation xsep(i) 
  % along the positive x-axis
  if strcmp(comsol,'3.5')
    fem = coli_model_35(xsep(i));
  else
    fem = coli_model_4x(xsep(i));
  end
  
  %run parameter sweep
  disp(['Running sequence ',num2str(i),'/',num2str(Nval)]);
  umod = urdme(fem,'huang', ...
                    'mode','bg',...
                    'verbose',1,...
                    'report',1,...
                    'outfile',sprintf('results/out%d.mat',i));

  % save input separately for later use
  if strcmp(comsol,'4.x')
    mphsave(umod.comsol,sprintf('results/mod%d.mph',i));
    umod=rmfield(umod,'comsol');
  end
  save(sprintf('results/in%d.mat',i),'umod');
end
