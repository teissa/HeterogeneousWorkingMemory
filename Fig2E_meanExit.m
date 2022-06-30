% papers (Bayes, Wimmer, etc)
clear all
close all

sigs=[0.01:0.01:0.1];

A=0.1;
ns=1:12;
figure; hold on;

clrs=[linspace(0,0.8,length(sigs));linspace(0,0.8,length(sigs));linspace(0,0.8,length(sigs))]';
for s=1:length(sigs)
%     D=sigs(s)^2/2
    sg=sigs(s);
Ttheor = (2*pi)^2*besseli(0,2*A/sg^2./ns)/sg^2./ns.^2;
    semilogy(ns,Ttheor,'linewidth',5,'color',clrs(s,:))
    drawnow; pause

%     for nn=1:length(ns)
%         n=ns(nn);
% %         In=(2*pi/n)*besseli(0,A/sg^2);
% %         MET(s,nn)=In^2/D;
% 
%     end
    %plot(ns,MET(s,:),'linewidth',5,'color',clrs(s,:))
end
 legend('0.01','0.02','0.03','0.04','0.05','0.06','0.07','0.08','0.09','0.1')   
 xlabel('Wells','Interpreter','Latex')
ylabel('Mean Exit Time','Interpreter','Latex')
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')

 %%
 clear;
 close;
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
