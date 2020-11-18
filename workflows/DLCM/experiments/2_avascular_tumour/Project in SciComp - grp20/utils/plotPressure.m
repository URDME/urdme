%% Plot the pressure in a 2D patched way
% Get the pressure for the active dofs
Pr_ = full(U); Pr_(adof) = Pr(adof_);

% Patch active dofs
ii = find(Pr_);
ax1 = axes;
patch(ax1, 'Faces',R(ii,:),'Vertices',V, ...
         'FaceVertexCData',Pr_(ii), 'FaceColor','flat');
hold on;
% Set colormap for the active dofs
map_start = graphics_color('bluish green');
map_stop = graphics_color('vermillion');
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([1.1*min(min(Pr_(ii))) 1.1*max(max(Pr_(ii)))]);
colormap(ax1, mymap)

% Get the pressure for the boundary dofs
Pr_(adof) = 0;
Pr_(idof) = Pr(idof_);

% Patch boundary dofs
ii = find(Pr_);
ax2 = axes;
patch(ax2, 'Faces',R(ii,:),'Vertices',V, ...
         'FaceVertexCData',Pr_(ii), 'FaceColor','flat');

% Set colormap for the boundary dofs
map_start = [0,0,1];
map_stop = [0.7,0.7,1];
xx = linspace(0,1,10);
map_matrix = map_start' + xx.*(map_stop' - map_start');
mymap = map_matrix';
caxis([1.1*min(min(Pr_(ii))) 1.1*max(max(Pr_(ii)))]);
colormap(ax2, mymap)

% Add quivers to the moving boundary and source dofs
% Compute all "moveb"-rates
x = zeros(1,0); y = zeros(1,0);
u = zeros(1,0); v = zeros(1,0);
k_ = 0;
for i = bdof_m'
  k_ = k_+1;
  jj_ = find(N(i,Adof));
  keep = find(U(Adof(jj_)) == 0);
  jj_ = jj_(keep);
  for j_ = jj_
    if Pr(bdof_m_(k_)) > Pr(j_)
      x = [x P(1,i)];
      y = [y P(2,i)];
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient*(Pr(bdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+1)*gradquotient*(Pr(bdof_m_(k_))-Pr(j_))];
    end
  end
end

% Compute all "moves"-rates
k_ = 0;
for i = sdof_m'
  k_ = k_+1;
  jj_ = find(N(i,Adof));
  keep = find(U(Adof(jj_)) < 2);
  jj_ = jj_(keep);
  for j_ = jj_
    if Pr(sdof_m_(k_)) > Pr(j_)
      x = [x P(1,i)];
      y = [y P(2,i)];
      u = [u (P(1,Adof(j_))-P(1,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
      v = [v (P(2,Adof(j_))-P(2,i))*Drate_(2*VU(Adof(j_))+U(Adof(j_))+1)*gradquotient*(Pr(sdof_m_(k_))-Pr(j_))];
    end
  end
end
ax3 = axes;
q = quiver(ax3, x,y,u,v);
q.Color = 'yellow';

%%Link them together and hide them
linkaxes([ax1,ax2,ax3])
axis equal;
for a = [ax1, ax2,ax3]
    a.Visible = 'off';
    a.XTick = [];
    a.YTick = [];
end
ax1.Visible = 'on';
%Then add colorbars and get everything lined up
cb1 = colorbar(ax1,'Position',[.09 .11 .02 .815]);
cb2 = colorbar(ax2,'Position',[.91 .11 .02 .815]);
title(ax1, sprintf('alpha = %d', alpha))