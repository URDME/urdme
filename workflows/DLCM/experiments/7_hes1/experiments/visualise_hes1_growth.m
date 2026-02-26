%VISUALISE_HES1_GROWTH Visualise Hes1-Notch dynamics
%   Plots snapshots of the dynamics for both continuous and discrete
%   models, using the data/hes1growth_... mat-files.

% E. Blom 2025-11-08

%% (1) continuous
% no proliferation
load hes1growth_cont_still.mat

mu_prols = [0, 1/(13*60), 1/(20*60)]; %which mu_prol values that are visualised

X = data;     % P
Smax = max(X(:));
X = X./Smax;  % normalize to [0,1]

figure()
t = tiledlayout(4,4, 'Padding', 'none', 'TileSpacing', 'none');
tidx = [1 33 50 100];

n = 1;      % tile counter
for i = tidx
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(X(:,i)>0);
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  if i == 1
    text(-0.6,0.65, "a) $\mu_{\mathrm{prol}}^{-1} = \infty$", ...
    'Interpreter','latex')
  else
    text(-0.3,0.65, "$t = " + i + "\%$", ...
    'Interpreter','latex')
  end
  clim([0 1])
  xmax = 1;
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

% half-fast proliferation
load hes1growth_cont_slow.mat

X = data;     % P
X = X./Smax;  % normalize to [0,1]

tidx = [1 33 50 100];

n = 1;        % tile counter
for i = tidx
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(X(:,i)>0);
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  if i == 1
    text(-0.6,0.65, "b) $\mu_{\mathrm{prol}}^{-1} = 20$h", ...
    'Interpreter','latex')
  end
  clim([0 1])
  xmax = 1;
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

% fast proliferation
load hes1growth_cont_fast.mat

X = data;     % P
X = X./Smax;  % normalize to [0,1]

tidx = [1 33 50 100];

n = 1;        % tile counter
for i = tidx
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(X(:,i)>0);
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  if i == 1
    text(-0.6,0.65, "c) $\mu_{\mathrm{prol}}^{-1} = 13$h", ...
    'Interpreter','latex')
  end
  clim([0 1])
  xmax = 1;
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

% then calculate the patterning...
nexttile([1 4])

load hes1growth_cont_samples.mat

samples = 10;

% create neighbor operators
N_op = dt_neighe(V,R);
N_op = N_op ~= 0;

% postprocess given save_P and save_P_uds
W_num = zeros(numel(mu_prol_all),1);
W_numCI = zeros(numel(mu_prol_all),2);
for i = 1:numel(mu_prol_all)
  % accumulator over NMC trials:
  Nsucc = 0;
  Ntrial = 0;
  for k = 1:samples
    P = reshape(mumod_U_all(i,k,:), [size(N_op,1),1]); % P = save_P(:,i,k);
    % assume that all with P = 0 are vacant voxels -- don't count these!
    Pdof = find(P>0);
    P = P(Pdof);
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

    % count once...
    neigh = N_op(Pdof,Pdof)*test;
    nhi_ = neigh(ixhi);
    % internal connections:
    NA = sum(nhi_(real(nhi_)+imag(nhi_) == 6));

    % ...count twice
    test(real(neigh)+imag(neigh) == 6) = 0;
    neigh = N_op(Pdof,Pdof)*test;
    nhi__ = neigh(ixhi);
    % boundary connections:
    NB = sum(nhi__(real(nhi_)+imag(nhi_) == 6));

    % successes/trials:
    nsucc = (imag(NA)-imag(NB))/2+imag(NB);
    ntrial = nsucc+real(NA);

    Nsucc = Nsucc+nsucc;
    Ntrial = Ntrial+ntrial;
    % done accumulating?
    if k == samples
      % then estimate:
      [W_num(i),W_numCI(i,:)] = binofit(Nsucc,Ntrial);
      if isnan(W_num(i))
        W_num(i) = 0; % assume its due to 0/0, which is 0 patterning
        W_numCI(i,:) = 0;
      end
    end
  end
end

% visualize
hold on
errorshade(mu_prol_all,W_numCI(:,1)',W_numCI(:,2)','r');
plot(mu_prol_all,W_num','LineWidth',2,'Color','r');
xline((mu_prols+1e-5),'k-.');  % hours^-1
yline(0.5,'k--');

xlabel(['proliferation rate $\mu_{\mathrm{prol}}$ ' ...
  '$\left[\mathrm{min}^{-1}\right]$'],'Interpreter','latex');
ylabel('patterning coefficient $p$','Interpreter','latex');
ylim([0.0 1]);
grid on

cb = colorbar;
set(cb,'Position',[0.034 0.38 .015 0.58], 'TickLabelInterpreter',...
        'latex');

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 420 400]);
set(gca, 'fontname', 'Roman', 'FontSize', 8.0)
set(gca,'TickLabelInterpreter',...
        'latex');

% tweak these parameters to fully remove figure whitespace...
t.Position = [0.053 0.055 0.931 0.95];% left, bottom, right, top

% uncomment to save:
%exportgraphics(t,'DN.pdf', 'resolution', 600)

%% (2) Plot discrete experiments
% no proliferation
load hes1growth_disc_still.mat

X = data;     % P
Smax = max(X(:));
X = X./Smax;  % normalize to [0,1]

figure()
t = tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'none');
tidx = [1 33 50 100];

n = 1;      % tile counter
for i = tidx
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(U(:,i)>0); % discrete levels can reach zero, so use U!
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  if i == 1
    text(-0.6,0.65, "a) $\mu_{\mathrm{prol}}^{-1} = \infty$", ...
    'Interpreter','latex')
  else
    text(-0.3,0.65, "$t = " + i + "\%$", ...
    'Interpreter','latex')
  end
  clim([0 1])
  xmax = 1;
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

% Proliferation = 1/(20*60)
load hes1growth_disc_slow.mat

X = data;     % P
X = X./Smax;  % normalize to [0,1]

tidx = [1 33 50 100];%[1 50 75 100];

n = 1;      % tile counter
for i = tidx
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  hold on
  Adof = find(U(:,i)>0); % discrete levels can reach zero, so use U!
  patch('Faces',R(Adof,:),'Vertices',V, 'FaceVertexCData', X(Adof,i), ...
  'FaceColor','flat', 'EdgeColor','none', 'Linewidth', 0.3);
  patch('Faces',R(Adof,:),'Vertices',V, ...
  'FaceColor','none', 'EdgeColor','k', 'Linewidth', 0.3);
  if i == 1
    text(-0.6,0.65, "b) $\mu_{\mathrm{prol}}^{-1} = 20$h", ...
    'Interpreter','latex')
  end
  clim([0 1])
  xmax = 1;
  axis(0.85*[-1 xmax -1 1])
  axis off
  n = n+1;
end

cb = colorbar;
set(cb,'Position',[0.034 0.1 .015 0.83], 'TickLabelInterpreter',...
        'latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 420 200]);
set(gca, 'fontname', 'Roman', 'FontSize', 8.0)
set(gca,'TickLabelInterpreter',...
        'latex');

% tweak these parameters to fully remove figure whitespace...
t.Position = [0.06 0.00 0.945 1.0]; % left, bottom, right, top
