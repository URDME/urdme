% Plot sum of V in different parts of the neuron domain. 

function xx = analyze_neuron(fem)

V  = 1;
Vk = 2;
Vd = 3;

Ms =  3;

totV = sum(fem.urdme.U);
%normalized time
tspan = fem.urdme.tspan/fem.urdme.tspan(end);

axon = fem.urdme.axon;
soma = fem.urdme.soma;
dendrites = fem.urdme.dendrites;

sV  = fem.urdme.U(1:Ms:end,:);
sVk = fem.urdme.U(2:Ms:end,:);
sVd = fem.urdme.U(3:Ms:end,:);

sA = sV+sVk+sVd;

max(sum(sA)./totV)

Vaxon      = sum(sA(axon,:));
Vsoma      = sum(sA(soma,:));
Vdendrites = sum(sA(dendrites,:));

figure('Position',[0 0 1200 600]);
set(gca,'FontSize',14)
plot(tspan,Vaxon./totV,tspan,Vsoma./totV,tspan,Vdendrites./totV,'-*','LineWidth',2);
% Plot vertical bar dividing the scenarios
hold
stem(0.5,1,'-black','Linewidth',2);
%annotation('textbox',[0.7 0.9 0.7 0.9],'String','\sigma_{dk}=10\sigma_{kd}');
xlabel('t (normalized)');
ylabel('Fraction of total V')
legend('Axon','Soma','Dendrites');



