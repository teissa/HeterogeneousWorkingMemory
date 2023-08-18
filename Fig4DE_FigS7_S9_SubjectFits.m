clear 
load Long_modelFits
subjects=120;
close all

Long_bestModel=BestModelMSE;
sum_longBest=histcounts(Long_bestModel);
LongMSE=MSE;


load Short_modelFits


Short_bestModel=BestModelMSE;
sum_shortBest=histcounts(Short_bestModel);
ShortMSE=MSE;


long_flat=sum(Long_bestModel==1); long_fixed=sum(Long_bestModel==2|Long_bestModel==3|Long_bestModel==4); 
long_learn=sum(Long_bestModel==5|Long_bestModel==6); long_learnDist=sum(Long_bestModel==7|Long_bestModel==8); 

short_flat=sum(Short_bestModel==1); short_fixed=sum(Short_bestModel==2|Short_bestModel==3|Short_bestModel==4); 
short_learn=sum(Short_bestModel==5|Short_bestModel==6); short_learnDist=sum(Short_bestModel==7|Short_bestModel==8); 

figure; hold on; bar([short_flat,short_fixed,short_learn,short_learnDist;long_flat,long_fixed,long_learn,long_learnDist],'stacked')
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xticks([1,2]); xticklabels({'Short','Long'})

figure; hold on; bar([sum_shortBest;sum_longBest],'stacked')
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xticks([1,2]); xticklabels({'Short','Long'})


Same_bestModel=Long_bestModel==Short_bestModel;

%model classes
Long_bestModelClass=Long_bestModel; 
Long_bestModelClass(Long_bestModelClass==2|Long_bestModelClass==3|Long_bestModelClass==4)=2;
Long_bestModelClass(Long_bestModelClass==5|Long_bestModelClass==6|Long_bestModelClass==7|Long_bestModelClass==8)=3;

Short_bestModelClass=Short_bestModel; 
Short_bestModelClass(Short_bestModelClass==2|Short_bestModelClass==3|Short_bestModelClass==4)=2;
Short_bestModelClass(Short_bestModelClass==5|Short_bestModelClass==6|Short_bestModelClass==7|Short_bestModelClass==8)=3;

Same_bestModelClass=Long_bestModelClass==Short_bestModelClass;
Diff_bestModelClass=Long_bestModelClass~=Short_bestModelClass;


%% offsets
offset_adjusted=offset; offset_adjusted(offset>90&offset<=180)=offset_adjusted(offset>90&offset<=180)-90;
offset_adjusted(offset>180&offset<=270)=offset_adjusted(offset>180&offset<=270)-180;offset_adjusted(offset>270&offset<=360)=offset_adjusted(offset>270&offset<=360)-270;

figure; hold on; histogram(offset_adjusted(Long_bestModel==2|Long_bestModel==3|Long_bestModel==4),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Short_bestModel==2|Short_bestModel==3|Short_bestModel==4),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Long_bestModel==5|Long_bestModel==6|Long_bestModel==7|Long_bestModel==8),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Short_bestModel==5|Short_bestModel==6|Short_bestModel==7|Short_bestModel==8),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Same_bestModelClass),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Diff_bestModelClass),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Long_bestModelClass(Same_bestModelClass)==3),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

figure; hold on; histogram(offset_adjusted(Long_bestModelClass(Same_bestModelClass)==2),0:5:90); plot([0 0],[0 20],'k--','LineWidth',5); ylim([0 20]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Subjects','Interpreter','latex'); xlabel('Offset','Interpreter','latex');

%histograms
twoDhist=histcounts2(Short_bestModelClass,Long_bestModelClass);

figure; hold on; imagesc(twoDhist); colormap("parula"); clim([0 30]); colorbar
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Short','Interpreter','latex'); xlabel('Long','Interpreter','latex'); axis('tight')
xticks([1,2,3]); xticklabels({'Flat','Fixed','Learn'})
yticks([1,2,3]); yticklabels({'Flat','Fixed','Learn'})