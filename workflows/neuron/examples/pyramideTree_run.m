%PYRAMIDETREE_RUN Complete coupling of eight neurons.
%   Starts a chain reaction by setting off a pulse in the "root"
%   neuron (#1). The full tree can be seen in
%   data/pyramide_tree.eps. An animated solution can be found in the
%   same folder as a .gif-animation.
%
%   The output is eight graphs, one for each neuron, showing the
%   membrane potential and synaptic current for that neuron over time.

% S. Engblom 2019-11-29 (Minor revision)
% S. Engblom 2018-06-21 (Revision)
% A. Senek 2017-05-31

if ~exist('report','var')
  report = 1;
end

load('data/pyramide_tree.mat');
prop_relpath_chnl = 'channel/';

% these will make up the nodes where the synaptic connections attach
% to the post-synaptic tree
input_nodes = [3 35 44 523 601 622 629 68];

% these will make up the nodes where the synaptic connections attach
% to the pre-synaptic tree
Vm_nodes = [1 30 42 521 599 620 627 66];

% we don't want the Crank-Nicholson to connect these nodes by the
% cable equation
tree.dA(2:3,:) = 0;
tree.dA(input_nodes,:) = 0;   
clear umod_chnl;
umod_chnl.private.model = tree;
nVoxels = tree.nVoxels;
nSyns = 8;
input_vec = zeros(1,nVoxels);
umod_chnl.private.I_inj = @(V,t,i)((t > 0).*input_vec)';

% setup neuron channels
t_setup = [0 0.01];
[~, ~, umod_chnl,] = neuron_solver(t_setup,umod_chnl,prop_relpath_chnl);

% setup synaptic channels
bogus_Vm = zeros(length(Vm_nodes),1)-75;
[umod_syn,~] = synaptic_solver([],t_setup,bogus_Vm,0);

% ODE stepping
dt = 0.05;
t_vec = 0:dt:15;
V_out = zeros(numel(t_vec)-1,nVoxels);
V_out(1,:) = umod_chnl.private.channels.Erest';
Iin = zeros(numel(t_vec),nSyns);

for i = 1:numel(t_vec)-1
  t_tmp = t_vec([i i+1]);
  % channel simulation
  [V_temp,~,umod_chnl] = neuron_solver(t_tmp,umod_chnl,prop_relpath_chnl);

  % synaptic simulation
  V_pre = V_temp(2,Vm_nodes)';
  [umod_syn,I_syn] = synaptic_solver(umod_syn,t_tmp,V_pre,0);

  % set the input current accordingly with a source term for t < 2 [ms]
  inj_vec = zeros(1,nVoxels);
  inj_vec(input_nodes) = I_syn(2,:);
  umod_chnl.private.I_inj = @(V,t,i)((t > 0).*inj_vec)'+ ...
      (1.4*(t < 2).*[0 0 1 zeros(1,nVoxels-3)]');
  I_tmp = umod_chnl.private.I_inj(1,t_tmp(1),1);
  Iin(i+1,:) = I_tmp(input_nodes);
  V_out(i+1,:) = V_temp(2,:);

  % report progress
  if mod(i,round(numel(t_vec)/10)) == 1 && report
    percent = floor(i/numel(t_vec)*100);
    str = ['[' num2str(percent) '%%]-'];
    fprintf(str);
  end
end
if(report), disp('[100%]'); end

%  output voltage plot
if exist('plotting_off','var') && plotting_off, return; end
if exist('yyaxis','file')
  figure 
  subplot(4,2,1)
  title('1')
  yyaxis right
  plot(t_vec,Iin(:,1))  
  hold on      
  yyaxis left
  plot(t_vec,V_out(:,576))
  
  subplot(4,2,2)
  title('4')
  yyaxis right
  plot(t_vec,Iin(:,4))  
  hold on           
  yyaxis left
  plot(t_vec,V_out(:,526))
  legend('Neuron potential','Input current');
  
  subplot(4,2,3)
  title('2')
  yyaxis right
  plot(t_vec,Iin(:,2))  
  hold on   
  yyaxis left
  plot(t_vec,V_out(:,438))
  
  subplot(4,2,4)
  title('5')
  yyaxis right
  plot(t_vec,Iin(:,5))  
  hold on 
  yyaxis left
  plot(t_vec,V_out(:,998))

  subplot(4,2,5)
  title('6')
  yyaxis right
  plot(t_vec,Iin(:,6))  
  hold on  
  yyaxis left
  plot(t_vec,V_out(:,828))

  subplot(4,2,6)
  title('7')
  yyaxis right
  plot(t_vec,Iin(:,7))  
  hold on    
  yyaxis left
  plot(t_vec,V_out(:,692)) 

  subplot(4,2,7)
  title('3')
  yyaxis right
  plot(t_vec,Iin(:,3))  
  hold on     
  yyaxis left
  plot(t_vec,V_out(:,363)) 

  subplot(4,2,8)
  title('8')
  yyaxis right
  plot(t_vec,Iin(:,8))  
  hold on   
  yyaxis left
  plot(t_vec,V_out(:,131))
else
  % use plotyy instead
  figure 
  subplot(4,2,1)
  title('1')
  plotyy(t_vec,Iin(:,1),t_vec,V_out(:,576))

  subplot(4,2,2)
  title('4')
  plotyy(t_vec,Iin(:,4),t_vec,V_out(:,526))
  legend('Neuron potential','Input current');
  
  subplot(4,2,3)
  title('2')
  plotyy(t_vec,Iin(:,2),t_vec,V_out(:,438))

  subplot(4,2,4)
  title('5')
  plotyy(t_vec,Iin(:,5),t_vec,V_out(:,998))

  subplot(4,2,5)
  title('6')
  plotyy(t_vec,Iin(:,6),t_vec,V_out(:,828))

  subplot(4,2,6)
  title('7')
  plotyy(t_vec,Iin(:,7),t_vec,V_out(:,692)) 

  subplot(4,2,7)
  title('3')
  plotyy(t_vec,Iin(:,3),t_vec,V_out(:,363)) 
  
  subplot(4,2,8)
  title('8')
  plotyy(t_vec,Iin(:,8),t_vec,V_out(:,131))
end

figure(1),
set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[300 400 680 550]); % for example

return;

% produce an animation
load('data/pyramide_tree.mat');
tree.R = tree.R';
make_gif(tree,t_vec,V_out,'pyramide_tree.gif');
