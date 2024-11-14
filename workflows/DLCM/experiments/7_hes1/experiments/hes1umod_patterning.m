%HES1UMOD_PATTERNING Pattern robustness for varying cell volume.
%
%   See also HES1UMOD, HES1UMOD2D_RUN.

% S. Engblom 2024-10-11 (Revision)
% G. Menz 2024-05-08

if ~exist('save_data','var')
  save_data = false;
end

% critical parameter to sweep over: volume of one cell/voxel
% vol has unit um^3 (so vol = 50 is ~one mouse stem cell)
vol = [1:4 5:5:50];
NMC = 10; % #independent replicas

% build the model using external script
solve = 0;
VOL = 1; % use unit scaling for the volume here
Nvoxels = 10;
hes1umod2D_run;

% modify the struct a bit:
umod.solve = 1; vmod.solve = umod.solve;
umod.parse = 1; vmod.parse = umod.parse;
umod.compile = 0; vmod.compile = umod.compile;
umod.report = 1; vmod.report = umod.report;
umod.tspan = linspace(0,36*60,10); vmod.tspan = umod.tspan;
rng(20240508);
umod.seed = randi(intmax('uint32'),NMC,1);

% output fate by the end
save_P = zeros(Nvoxels^2,numel(vol),NMC);

% (1) solve with deterministic interpretation first since it does not
% depend on the volume
disp('*** Pre-solving with UDS...');
vmod.u0 = u0*vol(1)/VOL;
vmod = urdme(vmod);
save_P_uds = vmod.U(5:Mspecies:end,end,:);
disp('   ...done.');

% (2) main solve: per volume and NMC i.i.d. replicas using NSM
for i = 1:numel(vol)
  disp(['*** Solving with volume = ' num2str(vol(i)) '...']);

  % scale the volume from previous round
  umod.vol = umod.vol*vol(i)/VOL;
  % note: use u0, not umod.u0 to avoid rounding errors:
  umod.u0 = repmat(round(u0*vol(i)/VOL),1,1,NMC);
  VOL = vol(i);

  % solve
  umod = urdme(umod);

  % postprocess and save end-result
  P = umod.U(5:Mspecies:end,end,:);

  % save end result
  save_P(:,i,:) = P;
  disp('   ...done.');
end

% postprocess given save_P and save_P_uds
W_num = zeros(numel(vol),1);
W_numCI = zeros(numel(vol),2);
for i = 0:numel(vol)
  for k = 1:NMC
    % accumulator over NMC trials:
    Nsucc = 0;
    Ntrial = 0;
    if i == 0
      % special compute for the uds-solver result
      P = save_P_uds;
    else
      P = save_P(:,i,k);
    end
    % find the largest increase in P
    [Pincr,ix] = sort(P);
    [~,ijmp] = max(diff(Pincr(1:floor(end*0.95))));
    % (except for possibly an end-effect at the 5% upper end)
    ixlo = ix(1:ijmp);
    ixhi = ix(ijmp+1:end);

    % count neighbors "lo-lo, lo-hi, hi-lo, hi-hi"
    test = zeros(size(P,1),1);
    test(ixlo) = 1; % "lo"
    test(ixhi) = 1i; % "hi" - imaginary unit to tell them apart
    neigh = N_op*test; % count neighbors as "lo+hi*1i"

    % restrict statistics to voxels with 6 neighbors
    neigh(real(neigh)+imag(neigh) ~= 6) = nan;

    % count once...
    neigh = N_op*test;
    nhi_ = neigh(ixhi);
    % internal connections:
    NA = sum(nhi_(real(nhi_)+imag(nhi_) == 6));

    % ...count twice
    test(real(neigh)+imag(neigh) == 6) = 0;
    neigh = N_op*test;
    nhi__ = neigh(ixhi);
    % boundary connections:
    NB = sum(nhi__(real(nhi_)+imag(nhi_) == 6));

    % successes/trials:
    nsucc = (imag(NA)-imag(NB))/2+imag(NB);
    ntrial = nsucc+real(NA);

    if i == 0
      % w_num_uds corresponds to coupling W(2,2), expected to be 1/2. We
      % estimate this in the form of a binomial fit, with #successes =
      % number of *unique* hi-hi couplings (by symmetricity of the
      % coupling relation, we get *twice* this number except for
      % boundary couplings), and #trials = number of potential
      % couplings (noting here that hi-lo couplings are uniquely
      % counted).
      [w_num_uds,w_num_udsCI] = binofit(nsucc,ntrial);
      break; % since NMC == 1 for uds
    else
      Nsucc = Nsucc+nsucc;
      Ntrial = Ntrial+ntrial;
      % done accumulating?
      if k == NMC
        % then estimate:
        [W_num(i),W_numCI(i,:)] = binofit(Nsucc,Ntrial);
      end
    end
  end
end

% data to save for being able to visualize
if save_data
  save('../data/urdme_patterning.mat','vol', ...
       'W_num','W_numCI','w_num_uds','w_num_udsCI', ...
       'save_P_uds','save_P');
end

% visualize
figure(1), clf, hold on,
% $$$ h = yline(w_num_uds,'b','LineWidth',2);
% $$$ col = get(h,'Color');
% $$$ errorshade([vol(1) vol(end)], ...
% $$$     w_num_udsCI([1 1]),w_num_udsCI([2 2]),col);
h = plot(vol,W_num','LineWidth',2,'Color','r');
col = get(h,'Color');
errorshade(vol,W_numCI(:,1)',W_numCI(:,2)',col);
plot(vol,W_num','LineWidth',2,'Color','r');
yline(0.5,'k--');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 280 180]);
set(gca,'TickLabelInterpreter',...
        'latex');
xlabel('Cell Volume $\left[\mu m^3\right]$','Interpreter','latex');
ylabel('$p$','Interpreter','latex');
ylim([0.3 0.7]);
box on

% make plot publishable:
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[200 200 260 160]);
%print -depsc ~/Desktop/rdme_patterning.eps
