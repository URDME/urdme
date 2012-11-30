% Run the annihilation example
%

umod.name = 'annihilation';
umod.test=1;

% Execute the solver. This will autimatically generate a propensity file
% from the inline proponsties described in the function annihilation.m
umod = urdme(umod,@annihilation);

Mspecies=2;
N = size(umod.mesh.p,2);
Z = umod.mesh.p(3,:);
%dofs=xmeshinfo(fem,'out','dofs');
%Z  = dofs.coords(3,1:Mspecies:end);

vol=umod.vol;

mean_A = mean(umod.U(1:2:end,50:101),2)./vol;
mean_B = mean(umod.U(2:2:end,50:101),2)./vol;

figure(1);plot(Z,mean_A,'-dr');hold on;
plot(Z,mean_B,'-db');hold off
legend('a','b','Location','N');
axis tight;