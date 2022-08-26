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
 edges=0:5:90;

%%
load LongTrialData
ofst=offset;
indof=find(ofst>90 &ofst<=180);
ofst(indof)=ofst(indof)-90;
indof=find(ofst>180 &ofst<=270);
ofst(indof)=ofst(indof)-180;
indof=find(ofst>270 &ofst<=360);
ofst(indof)=ofst(indof)-270;

 hL=histcounts(BestModelMSE);

  for k=2:3
     if k==2
        ind=find(BestModelMSE==2 |BestModelMSE==3 |BestModelMSE==4  );
        hHetL=histcounts(ofst(ind),edges);
     elseif k==3
        ind=find(BestModelMSE==5 |BestModelMSE==6 |BestModelMSE==7 |BestModelMSE==8 );
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

 hS=histcounts(BestModelMSE);
 
 for k=2:3
     if k==2
        ind=find(BestModelMSE==2 |BestModelMSE==3 |BestModelMSE==4 );
        hHetS=histcounts(ofst(ind),edges);
     elseif k==3
        ind=find(BestModelMSE==5 |BestModelMSE==6 |BestModelMSE==7 |BestModelMSE==8 );
        hLearnS=histcounts(ofst(ind),edges);
     end
 end
 
 %% summary plot
 
  figure; hold on; bar([hS;hL],'stacked')
 xticks([1 2]); xticklabels({'Short','Long'});
 ylabel('Subjects','Interpreter','Latex')
ylim([0 120])
xlim([0 4]); 
legend('Hom','Het 1','Het 2','Het Off','Learn Flat','Learn Het','Learn Dist Flat','Learn Dist Het','Interpreter','Latex')
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')


%% offsets

close all

figure; hold on
bar(edges(1:end-1),hLearnL);
 plot([0,0],[0,20],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 20])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

figure; hold on
bar(edges(1:end-1),hHetL);
 plot([0,0],[0,20],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 20])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')


figure; hold on
bar(edges(1:end-1),hLearnS);
 plot([0,0],[0,20],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 20])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')

figure; hold on
bar(edges(1:end-1),hHetS);
 plot([0,0],[0,20],'k--','linewidth',2)
xlabel('Offset','Interpreter','Latex'); ylabel('Subject Fraction','Interpreter','Latex')
ylim([0 20])
%xlim([0 nTrough1(end)]); 
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')