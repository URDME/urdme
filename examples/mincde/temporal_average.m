%  Postprocessing script for the mincde model. 
%  
%  Computes and plots the temporal average of MinD_m
%  along with the total number of MinD_m in one half of the
%  bacterium. 
%
%  INPUT: 
%      umod - An URMDE 1.2 model object with an attached solution.
%
%  P. Bauer,     2012-09-3
%  A. Hellander, 2010-06-9. 
%
%

function umod2 = temporal_average(umod)

pm = umod.sd == 2;
Mspecies =5;

% Temporal average of the populations
U = umod.U;
tspan    = umod.tspan;

U2 = sum(U,2)/(tspan(end)-tspan(1));

umod2=umod;
umod2.tspan = [];
U2 = [umod.U(:,1) U2];

umod2 = urdme2comsol(umod2,U2,[tspan(1) tspan(end)]);
if iscmp4x(umod.comsol) %Comsol 4.x
  xmi = mphxmeshinfo(umod.comsol);
  dofs = xmi.dofs; 
else
  dofs = xmeshinfo(umod.comsol,'Out','dofs');
end
 
[~,Ncells] = size(umod.comsol.mesh.p);

% 1D version
minz = min(umod.comsol.mesh.p(3,:));
maxz = max(umod.comsol.mesh.p(3,:));

vpm     = zeros(1,Ncells);
vpm(pm) = umod.vol(pm);
z  = dofs.coords(3,1:Mspecies:end);
zz = zeros(1,Ncells);
zz(pm) = z(pm); 

slices =linspace(minz,maxz,40);
MinDm = U2(2:Mspecies:end,end);

for i=1:numel(slices)
   ind = find(zz < slices(i));   
   expr = ['MinD_m*(z<' num2str(slices(i)) ')'];
   aver(i)=sum(MinDm(ind));
   vol(i)=sum(vpm(ind));      
end

MinD_c_atp  = umod.U(1:Mspecies:end,:);
MinD_c_adp  = umod.U(5:Mspecies:end,:);
MinE        = umod.U(3:Mspecies:end,:);
MinD_m      = umod.U(2:Mspecies:end,:);
MinDE       = umod.U(4:Mspecies:end,:);

ind1 = z <= (max(z)+min(z))/2;
ind2 = z >  (max(z)+min(z))/2;

tMinD_m1 = sum(MinD_m(ind1,:));
tMinD_m2 = sum(MinD_m(ind2,:));
tMinDE   = sum(MinDE(ind1,:)); 

tspan = umod.tspan;
meanspec = [mean(sum(MinD_c_atp)) mean(sum(MinD_m)) mean(sum(MinE)) mean(sum(MinDE)) mean(sum(MinD_c_adp))]
subplot(2,1,1); plot(tspan,tMinD_m1,'Linewidth',2);
title('Copy number in one of the poles. ');
xlabel('t [s]');
ylabel('# molecules');


conc = diff(aver)./diff(vol);
xline = 2.25e-6*ones(1,40);
vline = linspace(0,1,40); 
subplot(2,1,2); plot(slices(1:end-1)+0.5e-6,conc./max(conc),'*',xline,vline,'--','LineWidth',2);
set(gca,'YTick',[0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0],'YLim',[0 1],'XLim',[0 4.5e-6],'XTick',[0 2.25e-6 4.5e-6]);
title('Spatiotemporal average');

figure(2); postplot(umod2.comsol,'Tetdata','MinD_m');
title('Spatiotemporal average');
