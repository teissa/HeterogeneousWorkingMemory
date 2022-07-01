clear
close all

load LongTrialData


figure; histogram(BestModelKL);
xticks([1:4]); ylabel('Subject Count','interpreter','latex')
xticklabels({'hom','het','learn','learn dist'}); xlim([0 5]); ylim([0 120])
set(gca,'FontSize',30);set(gca,'TickLabelInterpreter','Latex')


    
flatKL=find(BestModelKL==1);
hetKL=find(BestModelKL==2);
learnKL=find(BestModelKL==3);
learndistKL=find(BestModelKL==4);
 hH=histcounts(offset(hetKL)*(pi/180),offEdges*(pi/180),'Normalization','probability');
 hL=histcounts(offset(learnKL)*(pi/180),offEdges*(pi/180),'Normalization','probability');
 hLD=histcounts(offset(learndistKL)*(pi/180),offEdges*(pi/180),'Normalization','probability');
[htest,pv]=ttest(hH,hLD);
figure; hold on; 
subplot(1,3,1);
xlabel('offset','interpreter','latex');set(gca,'FontSize',20);set(gca,'TickLabelInterpreter','Latex'); 
subplot(1,3,2);
xlabel('offset','interpreter','latex'); set(gca,'FontSize',20);set(gca,'TickLabelInterpreter','Latex'); 
subplot(1,3,3);
xlabel('offset','interpreter','latex');set(gca,'FontSize',20);set(gca,'TickLabelInterpreter','Latex');
