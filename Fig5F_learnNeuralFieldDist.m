clear
close all


%% plot total mean distortion for the neural field learning model and static models
load NeuralLearningCompare
% dlearn(2)=0.5;
% 
% for k=1:Nsim
%     mdLearn(k)=nanmean(dlearn(1:k));
% end
figure; hold on;
for k=1:10
    plot(mdLearn{k},'Color',[0.5,0.5,0.5],'LineWidth',2)
    mLearn(k,:)=mdLearn{k};
    mHom(k,:)=mdHom{k};
    mHet(k,:)=mdHet{k};
end
plot(nanmean(mLearn),'k','Linewidth',4)


plot(nanmean(mHom),'k--','linewidth',4);
plot(nanmean(mHom),'k--','linewidth',4);
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex')
xlabel('Trial','interpreter','latex','fontsize',30);
ylabel('Distortion','interpreter','latex','fontsize',30);


