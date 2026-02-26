%VISUALISE_CELL_SORTING Visualise cell sorting dynamics.
%   Plots snapshots from simulation data

% E. Blom 2024-11-19

%% Exp 1, 1st row:
load cellsort1.mat

U = umod.U;

tspan = [0 0.6 2.0 4.0]; %.*10^5;

tmap = 4200; % timescaling to hours
tframe = tspan*10^5/tmap; % time at snapshots below

zoomf = 0.75; % domain zoom-in factor
nrows = 3+1;
str = ['a)', 'b)', 'c)', 'd)'];
n = 1;  % loop counter
figure()
t = tiledlayout(nrows,4, 'Padding', 'none', 'TileSpacing', 'none');
for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  if i == 1 % write out cell type ratio
      text(-0.92-0.28,0,0, "$50\%$", ...
        'Interpreter','latex')
      text(-0.99-0.28,-0.12,0, "purple", ...
        'Interpreter','latex')
  end
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1+0.2 1+0.2])
  text(-0.7,+0.75,0, "$"+str(n:n+1)+"$", 'Interpreter','latex') % abcd)
  text(-0.35,+0.75,0, "$t="+round(tframe(round(n/2)))+"$ hours", ...
    'Interpreter','latex') % time
  axis off
  n = n+2;
end

%% Exp 2, 2nd row
load cellsort2.mat

U = umod.U;

for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  if i == 1 % write out cell type ratio
      text(-0.92-0.28,0,0, "$25\%$", ...
        'Interpreter','latex')
      text(-0.99-0.28,-0.12,0, "purple", ...
        'Interpreter','latex')
  end
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1 1])
  plot(-1:1, 0.65*[1 1 1], 'k--', 'linewidth', 1.0)
  axis off
end

%% Exp 3, 3nd row
load cellsort3.mat

U = umod.U;

for i = [1 15 50 100]
  nexttile
  patch('Faces',R(:,:),'Vertices',V, ... % grey background
  'FaceColor',[0.8 0.8 0.8], 'EdgeColor','none');
  if i == 1 % write out cell type ratio
      text(-0.92-0.28,0,0, "$13\%$", ...
        'Interpreter','latex')
      text(-0.99-0.28,-0.12,0, "purple", ...
        'Interpreter','latex')
  end
  hold on
  patch('Faces',R(U(1,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('sky blue'), 'EdgeColor','none');
  patch('Faces',R(U(2,:,i)>0,:),'Vertices',V, ...
  'FaceColor',graphics_color('reddish purple'), 'EdgeColor','none');
  axis(zoomf*[-1 1 -1-0.2 1-0.2])
  axis off
end

nexttile([1 4])
% calculate and plot fractional length
for n = 1:3
  load("cellsort"+n+".mat")
  % Get cell count from data
  U = reshape(sum(umod.U(1:2,:,:),1), ...
  size(umod.U,2),size(umod.U,3));
  % cell type data
  X = umod.U(:,:,:);

  % get Neighboring matrix from P, T
  nquants = 1; % number of field states pressure and nutrient
  Dexpr = cell(1,nquants+2);
  Dexpr(:) = {1};
  umod_tmp = pde2urdme(P,T,Dexpr);
  D = -umod_tmp.D(1:2+nquants:end, 1:2+nquants:end)';
  Ne = (D-diag(diag(D))~=0);

  [nani, nanj] = find(isnan(R)); % use to filter out NaNs

  T = 100; % = numel(umod.tspan)
  % calculate fractional length
  % total length of edges between voxels with different types:
  frac_len = zeros(T,1);
  n_edges = zeros(T,1); % keep track of full edge length each t
  for tt = 1:T
  adof = find(U(:,tt)>0);
  for v = adof' % adofs ordered the same as voxel idx in R
  ridx = 1:6;   % indices of R(v,:)
  if ismember(v,nani)
  ridx = 1:nanj-1;
  end
  ctype = find(X(:,v,tt)); % voxel celltype
  % find neighboring celltypes...
  [i, j] = find(X(:,find(Ne(v,:)),tt));
  het_frac = sum(i~=ctype)/6; % fraction of heterotypic voxel edge
  % we want the _total_ edge length of cells against different type
  tot_edge = perimeter(polyshape(V(R(v,ridx),1),V(R(v,ridx),2)));
  frac_len(tt) = frac_len(tt) + het_frac*tot_edge;
  % Count edges:
  % assuming equal length sides on each voxel and 6 edges per voxel (hex...):
  nactive = numel(i); % => # edges that will be counted _once_ more...
  % ... half the contribution of these edges, and neglect the edges to void
  % (no possibility for the latter to be heterotypic edges!):
  n_edges(tt) = n_edges(tt) + nactive/2;
  end
  frac_len(tt) = frac_len(tt)/2; % every edge between heterotypic voxels counted twice!
  % Note on occupied voxels: Edge length is still the same, since homotypic
  end
  hold on
  % assuming hex voxels here! (*6, etc))
  [frac, var] = binofit(round(frac_len/tot_edge*6),n_edges,0.01);
  plot(linspace(0, 4*10^5/tmap, 100), frac, 'linewidth', 2)
  errorshade(linspace(0, 4*10^5/tmap, 100), var(:,1)', var(:,2)')
  %plot(linspace(0, 4*10^5/tmap, 100), frac_len./max(frac_len), 'linewidth', 2)
end

grid on
% indicate when snapshots are taken
plot([tframe; tframe], [zeros(1,4); ones(1,4)], '--k', 'linewidth', 1.5)
ylabel('$\xi/\xi_0$', 'Interpreter', 'Latex')
xlabel('time [hours]', 'Interpreter', 'Latex')
axis([0 96 0 0.75])
legend('Case 1', 'Case 2', 'Case 3', 'Interpreter', 'Latex')

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[100 100 400 90*nrows]);
set(gca, 'fontname', 'Roman', 'FontSize', 10.0)
set(gca,'TickLabelInterpreter',...
        'latex');

% tweak these parameters to fully remove figure whitespace...
%t.Position = [+0.00 -0.075 1 1.08]; % left, bottom, right, top
% 2-by-4: t.Position = [0 -0.04 1 1.08];
% 3-by-4: t.Position = [0 -0.03 1 1.06];

% uncomment to save:
%exportgraphics(t,'cellsort1.pdf')