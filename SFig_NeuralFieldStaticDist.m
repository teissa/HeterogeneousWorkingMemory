%% Plot the mean distortion for the static neural field models

load NeuralLearningCompare

xboot=1000;
    for j=1:xboot
        r=randi(100,[1,100]);
        mdisthetBoot(j)=nanmean(dhet(r));
        mdisthomBoot(j)=nanmean(dhom(r));
        mdistlearnBoot(j)=nanmean(dlearn(r));
    end
    mdistHet=nanmean(mdisthetBoot);
    sdistHet=nanstd(mdisthetBoot);
    mdistHom=nanmean(mdisthomBoot);
    sdistHom=nanstd(mdisthomBoot);
    mdistLearn=nanmean(mdistlearnBoot);
    sdistLearn=nanstd(mdistlearnBoot);
figure; hold on;
bar(1,mdistHom,'k')
errorbar(1,mdistHom,sdistHom,'k.')
bar(2,mdistHet,'k')
errorbar(2,mdistHet,sdistHet,'k.')
% bar(3,mdistHet,'k')
% errorbar(3,mdistHet,sdistHet,'k.')
ylabel('$\bar{d}_{\rm{Tot}}$','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
 xlim([0,3]);xticks([1,2]); xticklabels([{'Hom','Het'}])
ylim([0 0.2])
