clear
close all

load hetP_varDT

D=0.05;
figure; hold on; 
clrs=[linspace(0,0.8,4);linspace(0,0.8,4);linspace(0,0.8,4)]';
for k=1:3
    T=Ts(k);
    
    % find distortion with flat poential (no troughs)

    xfs = linspace(-pi,pi,201); %angles discretized to find distortion
    dx = xfs(2)-xfs(1); %angle step size

    fsol = 1/2/pi; %baseline value

    for n=1:10 %solve using first 10 terms of Fourier series solution sum_{n=1}^{\infty} cos(n*x)*exp(-D*n^2*t)/pi+1/2/pi
        fsol = fsol+cos(n*xfs).*exp(-D*n^2*T)/pi; 
    end
    fdist = dx*sum((1-cos(xfs)).*fsol); %distortion compares the distances from the true value and sums them up, noramlized by descritization step size

    plot(rdist{2,k}(1:16,11)-fdist,'linewidth',5,'color',clrs(k,:))
end
plot([0 16],[0 0],'k--','linewidth',5)
plot([4 4],[-1 1],'r','linewidth',1)

set(gca,'fontsize',30);
set(gca,'ticklabelinterpreter','latex')
xlabel('Wells','interpreter','latex');
ylabel('Relative Distortion','interpreter','latex')
legend('1s','5s','10s','interpreter','latex')
ylim([-1 1]);
xlim([0 16])