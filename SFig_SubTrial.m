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

Same_bestModel=Long_bestModel==Short_bestModel;

% het model
D_H= 0.01:0.01:0.2;
nTrough=1:12;
A=0.1:0.1:2; %amplitude of attractors
N=1000; %number of spatial steps


%learning models
D_L = 0.01:0.01:0.2;   % noise corresponding to corrosion of parametric WM rep, diffusion
k_P =1:1:10;  % amplitude of von mises to be used
s_H =1/4;   % shift
s_C = 1:1:10;     % scaling

xfs = linspace(-pi,pi,360); dx = xfs(2)-xfs(1);
Ntrials=200;
dt = 0.01;  % timestep of stochastic sims

n=4;
A=1;
ut=@(x)-A*cos(n*x);
baselinePot_flat=1/2/pi+0*xfs; %initial potential landscape (flat)
baselinePot_het=1/2/pi+ut(xfs); %initial potential landscape (flat)

sameModels=Long_bestModel(Same_bestModel);
indSameModels=find(Same_bestModel==1);
indSameLearnModels=indSameModels(sameModels==5|sameModels==6|sameModels==7|sameModels==8);

%get MSEs and index for long trials
load Long_modelFits

for k=1:length(indSameLearnModels)
    [~,mt(k)]=min(LongMSE(indSameLearnModels(k),:));
    if mt(k)==5
        learnMSEcv(k,1:5)=lfMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,1:5)=lfMSEind(:,indSameLearnModels(k));
    elseif mt(k)==6
        learnMSEcv(k,1:5)=lhMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,1:5)=lhMSEind(:,indSameLearnModels(k));
    elseif mt(k)==7
        learnMSEcv(k,1:5)=ldfMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,1:5)=ldfMSEind(:,indSameLearnModels(k));
    elseif mt(k)==8
        learnMSEcv(k,1:5)=ldhMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,1:5)=ldhMSEind(:,indSameLearnModels(k));
    end

    hetMSEcv(k,1:5)=h1MSE(:,indSameLearnModels(k));
    hetMSEcv_ind(k,1:5)=h1MSEind(:,indSameLearnModels(k));
end

%get MSEs and index for short trials
load Short_modelFits

for k=1:length(indSameLearnModels)
    [~,mt(k)]=min(ShortMSE(indSameLearnModels(k),:));
    if mt(k)==5
        learnMSEcv(k,6:10)=lfMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,6:10)=lfMSEind(:,indSameLearnModels(k));
    elseif mt(k)==6
        learnMSEcv(k,6:10)=lhMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,6:10)=lhMSEind(:,indSameLearnModels(k));
    elseif mt(k)==7
        learnMSEcv(k,6:10)=ldfMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,6:10)=ldfMSEind(:,indSameLearnModels(k));
    elseif mt(k)==8
        learnMSEcv(k,6:10)=ldhMSE(:,indSameLearnModels(k));
        learnMSEcv_ind(k,6:10)=ldhMSEind(:,indSameLearnModels(k));
    end

    hetMSEcv(k,6:10)=h1MSE(:,indSameLearnModels(k));
    hetMSEcv_ind(k,6:10)=h1MSEind(:,indSameLearnModels(k));
end

figure(1); hold on; plot([0 200],[0 0],'k')
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('LLR','Interpreter','latex'); xlabel('Trial','Interpreter','latex');

% find best model overall
for k=1:length(indSameLearnModels)
    ind=(indSameLearnModels(k)-1)*Ntrials+1:indSameLearnModels(k)*Ntrials;
    Tar=wrapToPi(Target(ind) .*(pi/180));
    Dis=wrapToPi(Distractor(ind) .*(pi/180));
    off=wrapToPi(offset(k).*(pi/180));
    Tim=Time(ind)/1000;


    cResp=wrapToPi(Report(ind).*(pi/180));
    [~,mtL]=min(learnMSEcv(k,:));
    paramL=learnMSEcv_ind(k,mtL);
    %build learning potentials
    [kp,sc,d]=ind2sub([length(k_P),length(s_C),length(D_L)],paramL);
    
    if mt(k)==5
        learnPot = LearnParticleLandscape(Ntrials,baselinePot_flat,dx,xfs,Tar, s_C(sc), s_H, k_P(kp));
    elseif mt(k)==6
        learnPot = LearnParticleLandscape(Ntrials,baselinePot_het,dx,xfs,Tar, s_C(sc), s_H, k_P(kp));
    elseif mt(k)==7
        learnPot = LearnDistParticleLandscape(Ntrials,baselinePot_flat,dx,xfs,Tar, s_C(sc), s_H, k_P(kp), Dis);
    elseif mt(k)==8
        learnPot = LearnDistParticleLandscape(Ntrials,baselinePot_het,dx,xfs,Tar, s_C(sc), s_H, k_P(kp), Dis);
    end

     figure; imagesc(learnPot'); colorbar
     set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('$\theta$','Interpreter','latex'); xlabel('Trial','Interpreter','latex'); axis('tight')

  
    %%%% het potential 1 test
    [~,mtH]=min(hetMSEcv(k,:));
    paramH=hetMSEcv_ind(k,mtH);
    [dM,nnM,aaM]=ind2sub([length(D_H),length(nTrough),length(A)],paramH);

      %%%% pde 
    [bpostL(k,:),msesL(k,:),probL] = learnPDEsim_be(D_L(dM),learnPot, dt,Tim, N,Tar,Ntrials,cResp);

    [bpostH(k,:),msesH(k,:),probH] = hetPDEsim_be(D_H(dM),A(aaM),nTrough(nnM),dt,Tim, N,Tar,Ntrials,cResp,off);

    LLR(k,:)=log(bpostL(k,:)./bpostH(k,:));
    figure(1); plot(LLR(k,:),'LineWidth',2); drawnow;
end

for k=1:Ntrials
    fracLM(k)=sum(LLR(:,k)>0)/length(indSameLearnModels);
    totalLLRfracLM(k)=sum(runsumLLR(:,k)>0)/length(indSameLearnModels);
end
figure; hold on; plot(fracLM); ylim([0 1]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Fraction Learning Subjects','Interpreter','latex'); xlabel('Trial','Interpreter','latex'); 

figure; hold on; plot(totalLLRfracLM); ylim([0 1]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Fraction Learning Subjects','Interpreter','latex'); xlabel('Trial','Interpreter','latex'); 


figure; hold on;
for k=1:length(indSameLearnModels)
    plot(runsumLLR(k,:),'LineWidth',2);
end
plot([0 Ntrials],[0 0],'k--','LineWidth',5); %ylim([0 1]); 
set(gca,'FontSize',24); set(gca, 'TickLabelInterpreter','latex'); ylabel('Total LLR','Interpreter','latex'); xlabel('Trial','Interpreter','latex'); 
