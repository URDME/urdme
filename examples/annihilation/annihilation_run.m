% Run the annihilation example
%
clc;

umod.name = 'annihilation';
%umod.test=1;

% Execute the solver. This will autimatically generate a propensity file
% from the inline proponsties described in the function annihilation.m
%umod = urdme(umod,@annihilation);
%umod = urdme(umod,'annihilation');
%umod = urdme(umod,'annihilation',{'Solver','dfsp','verbose',1});
%umod = urdme(umod,'annihilation',{'Solver','dfsp','tau',1e-4,'max_jump',10,'verbose',0});
umod = urdme(umod,'annihilation',{'Solver','dfsp','verbose',1,'DFSP_cache','ann_dfsp_cache_file.mat'});


%setenv('URDME_SOLVER_PATH',strcat(strcat(pwd,'/Users/brian/Desktop/research/operator_splitting/urdme_solvers/:'),getenv('URDME_ROOT')))
%umod = urdme(umod,'annihilation',{'Solver','adaptive_dfsp'});

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