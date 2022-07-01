clear
close all

 %learning models
D_L = 0.01:0.01:0.2;   % noise corresponding to corrosion of parametric WM rep, diffusion
k_P =1:1:10;  % amplitude of von mises to be used
s_H =1/4;   % shift
s_C = 1:1:10;     % scaling
nParamL=3;
nParamLD=3;

% het model
D_H= 0.01:0.01:0.2;
nTrough1=1:12;
A1=0.1:0.1:2; %amplitude of attractors
nTrough2=1:12;
A2=0.1:0.1:2; %amplitude of attractors
offs=0:5:90; %offsets
offs=wrapToPi(offs.*(pi/180));
nParamH1=3;
nParamH2=5;
nParamHOff=4;


%flat model
D_F= 0.01:0.01:0.2;
nParamF=1;

%%
close all
load LongTrialData
ofst=offset;
indof=find(ofst>90 &ofst<=180);
ofst(indof)=ofst(indof)-90;
indof=find(ofst>180 &ofst<=270);
ofst(indof)=ofst(indof)-180;
indof=find(ofst>270 &ofst<=360);
ofst(indof)=ofst(indof)-270;

 hL=histcounts(BestModelMSE);
 figure; hold on;
 edges=0:5:90;
 for k=1:8
    ind=find(BestModelMSE==k);
    subplot(2,4,k); hold on;
    histogram(ofst(ind),edges)
    plot([0,0],[0,20],'r','linewidth',5)
%     plot([90,90],[0,1],'r','linewidth',5)
%     plot([180,180],[0,1],'r','linewidth',5)
%     plot([270,270],[0,1],'r','linewidth',5)
%     plot([360,360],[0,1],'r','linewidth',5)

    title(['model' num2str(k) 'Long'])
 end
 
 figure; hold on;
 for k=1:3
     if k==1
        ind=find(BestModelMSE==1);
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['hom long'])
     elseif k==2
        ind=find(BestModelMSE==2 |BestModelMSE==3 );
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['het long'])
        hHetL=histcounts(ofst(ind),edges);
     elseif k==3
        ind=find(BestModelMSE==5 |BestModelMSE==6 |BestModelMSE==7);
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['learn long'])
        hLearnL=histcounts(ofst(ind),edges);
     end
 end
        
 
 
 load ShortTrialData
ofst=offset;
indof=find(ofst>90 &ofst<=180);
ofst(indof)=ofst(indof)-90;
indof=find(ofst>180 &ofst<=270);
ofst(indof)=ofst(indof)-180;
indof=find(ofst>270 &ofst<=360);
ofst(indof)=ofst(indof)-270;

 hL=histcounts(BestModelMSE);
 figure; hold on;
 edges=0:5:90;
 for k=1:8
    ind=find(BestModelMSE==k);
    subplot(2,4,k); hold on;
    histogram(ofst(ind),edges)
    plot([0,0],[0,20],'r','linewidth',5)
%     plot([90,90],[0,1],'r','linewidth',5)
%     plot([180,180],[0,1],'r','linewidth',5)
%     plot([270,270],[0,1],'r','linewidth',5)
%     plot([360,360],[0,1],'r','linewidth',5)

    title(['model' num2str(k) 'short'])
 end
 
 figure; hold on;
 for k=1:3
     if k==1
        ind=find(BestModelMSE==1);
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['hom short'])
     elseif k==2
        ind=find(BestModelMSE==2 |BestModelMSE==3 );
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['het short'])
        hHetS=histcounts(ofst(ind),edges);
     elseif k==3
        ind=find(BestModelMSE==5 |BestModelMSE==6 |BestModelMSE==7);
        subplot(1,3,k); hold on;
    histogram(ofst(ind),edges)
        plot([0,0],[0,20],'r','linewidth',5)
%         plot([90,90],[0,1],'r','linewidth',5)
%         plot([180,180],[0,1],'r','linewidth',5)
%         plot([270,270],[0,1],'r','linewidth',5)
%         plot([360,360],[0,1],'r','linewidth',5)
        title(['learn short'])
        hLearnS=histcounts(ofst(ind),edges);
     end
 end

 
 
 %%
 load ShortTrialData
 hS=histcounts(BestModelMSE);
 figure; hold on;
 edges=0:20:360;
 for k=1:8
    ind=find(BestModelMSE==k);
    subplot(2,4,k); hold on;
    histogram(offset(ind),edges,'Normalization','probability')
    plot([0,0],[0,1],'r','linewidth',5)
    plot([90,90],[0,1],'r','linewidth',5)
    plot([180,180],[0,1],'r','linewidth',5)
    plot([270,270],[0,1],'r','linewidth',5)
    plot([360,360],[0,1],'r','linewidth',5)

    title(['model' num2str(k) 'Short'])
 end
 
 figure; hold on;
 for k=1:3
     if k==1
        ind=find(BestModelMSE==1);
        subplot(1,3,k); hold on;
        histogram(offset(ind),edges,'Normalization','probability')
        plot([0,0],[0,1],'r','linewidth',5)
        plot([90,90],[0,1],'r','linewidth',5)
        plot([180,180],[0,1],'r','linewidth',5)
        plot([270,270],[0,1],'r','linewidth',5)
        plot([360,360],[0,1],'r','linewidth',5)
        title(['hom short'])
     elseif k==2
        ind=find(BestModelMSE==2 |BestModelMSE==3);
        subplot(1,3,k); hold on;
        histogram(offset(ind),edges,'Normalization','probability')
        plot([0,0],[0,1],'r','linewidth',5)
        plot([90,90],[0,1],'r','linewidth',5)
        plot([180,180],[0,1],'r','linewidth',5)
        plot([270,270],[0,1],'r','linewidth',5)
        plot([360,360],[0,1],'r','linewidth',5)
        title(['het short'])
        hHetS=histcounts(offset(ind),edges,'Normalization','probability');
     elseif k==3
        ind=find(BestModelMSE==5 |BestModelMSE==6 |BestModelMSE==7);
        subplot(1,3,k); hold on;
        histogram(offset(ind),edges,'Normalization','probability')
        plot([0,0],[0,1],'r','linewidth',5)
        plot([90,90],[0,1],'r','linewidth',5)
        plot([180,180],[0,1],'r','linewidth',5)
        plot([270,270],[0,1],'r','linewidth',5)
        plot([360,360],[0,1],'r','linewidth',5)
        title(['learn short'])
        hLearnS=histcounts(offset(ind),edges,'Normalization','probability');
     end
 end
 

 
 figure; hold on; bar([hS;hL],'stacked')
 xticks([1 2]); xticklabels({'Short','Long'});
 ylabel('Subjects','Interpreter','Latex')
ylim([0 120])
xlim([0 4]); 
legend('Hom','Het 1','Het 2','Het Off','Learn Flat','Learn Het','Learn Dist Flat','Learn Dist Het','Interpreter','Latex')
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

        

 
 %%
 load LongTrialData
for k=1:8
    ind=find(BestModelMSE==k);
    
     if ~isempty(ind)
         for j=1:length(ind)
             if k==1
                 FminInd(j)=fMSEmin_ind{ind(j)};
                 dF(j)=D_F(FminInd(j));
             elseif k==3
                 H2minInd(j)=h2MSEmin_ind{ind(j)};
                 [dH2(j),nnH2_1(j),aaH2_1(j),nnH2_2(j),aaH2_2(j)]=ind2sub([length(D_H),length(nTrough1),length(A1),length(nTrough2),length(A2)],H2minInd(j));
                 dH2(j)=D_H(dH2(j));nnH2_1(j)=nTrough1(nnH2_1(j)); nnH2_2(j)=nTrough2(nnH2_2(j)); aaH2_1(j)=A1(aaH2_1(j)); aaH2_2(j)=A2(aaH2_2(j));  
             elseif k==4
                 HOffminInd(j)=hOffMSEmin_ind{ind(j)};
                 [dHO(j),nnHO(j),aaHO(j),ooHO(j)]=ind2sub([length(D_H),length(nTrough1),length(A1),length(offs)],HOffminInd(j));
                  dHO(j)=D_H(dHO(j));nnHO(j)=nTrough1(nnHO(j));aaHO(j)=A1(aaHO(j));ooHO(j)=offs(ooHO(j));
             elseif k==2
                 H1minInd(j)=h1MSEmin_ind{ind(j)};
                 [dH1(j),nnH1(j),aaH1(j)]=ind2sub([length(D_H),length(nTrough1),length(A1)],H1minInd(j));
                  dH1(j)=D_H(dH1(j));nnH1(j)=nTrough1(nnH1(j));aaH1(j)=A1(aaH1(j));
             elseif k==5
                 LFminInd(j)=lMSEfmin_ind{ind(j)};
                 [kpLF(j),scLF(j),dLF(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LFminInd(j));
                kpLF(j)=k_P(kpLF(j)); scLF(j)=s_C(scLF(j));dLF(j)=D_L(dLF(j));
             elseif k==6
                 LHminInd(j)=lMSEhmin_ind{ind(j)};
                 [kpLH(j),scLH(j),dLH(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LHminInd(j));
                kpLH(j)=k_P(kpLH(j)); scLH(j)=s_C(scLH(j));dLH(j)=D_L(dLH(j));
             elseif k==7
                 LdFminInd(j)=ldMSEfmin_ind{ind(j)};
                  [kpLdF(j),scLdF(j),dLdF(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LdFminInd(j));
                kpLdF(j)=k_P(kpLdF(j)); scLdF(j)=s_C(scLdF(j));dLdF(j)=D_L(dLdF(j));
             elseif k==8
                 LdHminInd(j)=ldMSEhmin_ind{ind(j)};
                 [kpLdH(j),scLdH(j),dLdH(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LdHminInd(j));
                kpLdH(j)=k_P(kpLdH(j)); scLdH(j)=s_C(scLdH(j));dLdH(j)=D_L(dLdH(j));

             end
         end
     end
end

%homogeneous
figure; hold on; histogram(dF,D_F,'Normalization','probability');
xlabel('Diffusion','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 1])
xlim([0 D_F(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

%heterogeneous
nTr=[nnH1,nnH2_1,nnH2_2,nnHO];
figure; hold on; histogram(nTr,nTrough1,'Normalization','probability');
xlabel('Wells','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.6])
xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')


%learning
KpSc=[kpLF.*scLF, kpLH.*scLH, kpLdF.*scLdF, kpLdH.*scLdH];
figure; hold on; histogram(KpSc,0:10:100,'Normalization','probability');
xlabel('Learning Parameter','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.6])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')


 %%
 load ShortTrialData
for k=1:8
    ind=find(BestModelMSE==k);
    
     if ~isempty(ind)
         for j=1:length(ind)
             if k==1
                 FminInd(j)=fMSEmin_ind{ind(j)};
                 dF(j)=D_F(FminInd(j));
             elseif k==3
                 H2minInd(j)=h2MSEmin_ind{ind(j)};
                 [dH2(j),nnH2_1(j),aaH2_1(j),nnH2_2(j),aaH2_2(j)]=ind2sub([length(D_H),length(nTrough1),length(A1),length(nTrough2),length(A2)],H2minInd(j));
                 dH2(j)=D_H(dH2(j));nnH2_1(j)=nTrough1(nnH2_1(j)); nnH2_2(j)=nTrough2(nnH2_2(j)); aaH2_1(j)=A1(aaH2_1(j)); aaH2_2(j)=A2(aaH2_2(j));  
             elseif k==4
                 HOffminInd(j)=hOffMSEmin_ind{ind(j)};
                 [dHO(j),nnHO(j),aaHO(j),ooHO(j)]=ind2sub([length(D_H),length(nTrough1),length(A1),length(offs)],HOffminInd(j));
                  dHO(j)=D_H(dHO(j));nnHO(j)=nTrough1(nnHO(j));aaHO(j)=A1(aaHO(j));ooHO(j)=offs(ooHO(j));
             elseif k==2
                 H1minInd(j)=h1MSEmin_ind{ind(j)};
                 [dH1(j),nnH1(j),aaH1(j)]=ind2sub([length(D_H),length(nTrough1),length(A1)],H1minInd(j));
                  dH1(j)=D_H(dH1(j));nnH1(j)=nTrough1(nnH1(j));aaH1(j)=A1(aaH1(j));
             elseif k==5
                 LFminInd(j)=lMSEfmin_ind{ind(j)};
                 [kpLF(j),scLF(j),dLF(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LFminInd(j));
                kpLF(j)=k_P(kpLF(j)); scLF(j)=s_C(scLF(j));dLF(j)=D_L(dLF(j));
             elseif k==6
                 LHminInd(j)=lMSEhmin_ind{ind(j)};
                 [kpLH(j),scLH(j),dLH(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LHminInd(j));
                kpLH(j)=k_P(kpLH(j)); scLH(j)=s_C(scLH(j));dLH(j)=D_L(dLH(j));
             elseif k==7
                 LdFminInd(j)=ldMSEfmin_ind{ind(j)};
                  [kpLdF(j),scLdF(j),dLdF(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LdFminInd(j));
                kpLdF(j)=k_P(kpLdF(j)); scLdF(j)=s_C(scLdF(j));dLdF(j)=D_L(dLdF(j));
             elseif k==8
                 LdHminInd(j)=ldMSEhmin_ind{ind(j)};
                 [kpLdH(j),scLdH(j),dLdH(j)]=ind2sub([length(k_P),length(s_C),length(D_L)],LdHminInd(j));
                kpLdH(j)=k_P(kpLdH(j)); scLdH(j)=s_C(scLdH(j));dLdH(j)=D_L(dLdH(j));

             end
         end
     end
end

%homogeneous
figure; hold on; histogram(dF,D_F,'Normalization','probability');
xlabel('Diffusion','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 1])
xlim([0 D_F(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

%heterogeneous
nTr=[nnH1,nnH2_1,nnH2_2,nnHO];
figure; hold on; histogram(nTr,nTrough1,'Normalization','probability');
xlabel('Wells','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.6])
xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')


%learning
KpSc=[kpLF.*scLF, kpLH.*scLH, kpLdF.*scLdF, kpLdH.*scLdH];
figure; hold on; histogram(KpSc,0:10:100,'Normalization','probability');
xlabel('Learning Parameter','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.6])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

%%


figure; hold on
bar(edges(1:end-1),hLearnL); 
 plot([0,0],[0,1],'k--','linewidth',2)
plot([90,90],[0,1],'k--','linewidth',2)
plot([180,180],[0,1],'k--','linewidth',2)
plot([270,270],[0,1],'k--','linewidth',2)
plot([360,360],[0,1],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.2])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

figure; hold on
bar(edges(1:end-1),hLearnS); 
 plot([0,0],[0,1],'k--','linewidth',2)
plot([90,90],[0,1],'k--','linewidth',2)
plot([180,180],[0,1],'k--','linewidth',2)
plot([270,270],[0,1],'k--','linewidth',2)
plot([360,360],[0,1],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 0.2])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')
