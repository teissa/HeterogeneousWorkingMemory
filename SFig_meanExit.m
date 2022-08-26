%Mean Exit Time

clear 
close all

sigs=[0.02:0.02:0.1];
clrs=[linspace(0,0.8,length(sigs));linspace(0,0.8,length(sigs));linspace(0,0.8,length(sigs))]';

for k=1:length(sigs)
 s = sigs(k);    % sigma (std) of the noise term in the SDE
A = 0.1;    % unscaled height of the periodic potential
ns = [1:12];    % frequency of the periodic potential

% plot the theory curve
Ttheor = (2*pi)^2*besseli(0,2*A/s^2./ns)/s^2./ns.^2;
semilogy(ns,Ttheor,'Color',clrs(k,:),'linewidth',5);
hold on
set(gca,'fontsize',24)
xlabel('$n$','fontsize',30,'interpreter','latex')
ylabel('$T$','fontsize',30,'interpreter','latex')

end
 legend('0.02','0.04','0.06','0.08','0.1')   
 xlabel('Wells','Interpreter','Latex')
ylabel('Mean Exit Time','Interpreter','Latex')
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
