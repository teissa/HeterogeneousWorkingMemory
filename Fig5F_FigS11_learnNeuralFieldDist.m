clear
close all


%% plot total mean distortion for the neural field learning model and static models
load GaussLearnCE_decay

figure; hold on;
plot(mdLearn,'Color',[0.5,0.5,0.5],'LineWidth',3)
plot([0 Nsim],[mdHom(end), mdHom(end)],'k--','LineWidth',3)
plot([0 Nsim],[mdHet(end), mdHet(end)],'k.-','LineWidth',3)
set(gca,'FontSize',24)
xlabel('trial')
ylabel('mean distortion')

%% changing environment

figure; hold on;
plot(mdLearnCE,'Color',[0.5,0.5,0.5],'LineWidth',3)
plot([0 Nsim],[mdHom(end), mdHom(end)],'k--','LineWidth',3)
plot([0 Nsim],[mdHet(end), mdHet(end)],'k.-','LineWidth',3)
set(gca,'FontSize',24)
xlabel('trial')
ylabel('mean distortion')

figure; hold on; 
pcolor(xsc,xsc,W_learnCE); shading flat,  axis('tight')
title('Ending Weights','Interpreter','Latex'); plot([0 0],[-180 180],'k'); plot([90 90],[-180 180],'k');plot([-90 -90],[-180 180],'k')
caxis([-0.01 0.05]); colorbar
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex')
xlabel('Preferred Target (degrees)','Interpreter', 'latex')
ylabel('Preferred Target (degrees)','Interpreter', 'latex')



%starting weights
hetCE=cos(pi*n*Y/180+shft);
W_learnCE= dx*(1+s*hetCE).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));


figure; hold on; 
pcolor(xsc,xsc,W_learnCE); shading flat,  axis('tight')
title('Starting Weights','Interpreter','Latex'); plot([0 0],[-180 180],'k'); plot([90 90],[-180 180],'k');plot([-90 -90],[-180 180],'k')
caxis([-0.01 0.05]); colorbar
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex')
xlabel('Preferred Target (degrees)','Interpreter', 'latex')
ylabel('Preferred Target (degrees)','Interpreter', 'latex')