clear
close all
load HumDataExp2
%% parameter distributions

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

%
nSim=100;
nParam=100;
edges=linspace(-pi,pi,36);
xfs = linspace(-pi,pi,360); dx = xfs(2)-xfs(1);
Trials=100;
T = 4;      % delay time (seconds)
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim

n=4;
A=1;
ut=@(x)-A*cos(n*x);
baselinePot_flat=1/2/pi+0*xfs; %initial potential landscape (flat)
baselinePot_het=1/2/pi+ut(xfs); %initial potential landscape (flat)
lambdaMSE=0.1; %penalty scaling for complexity in MSE
lambdaKL=0.01; %penalty scaling for complexity in KL
%% testing by subject for short mixed trials
subjects=120;cv=1;
fMSEmin_ind=cell(1,subjects);fKLmin_ind=fMSEmin_ind;
h1MSEmin_ind=fMSEmin_ind; h1KLmin_ind=fKLmin_ind;
h2MSEmin_ind=fMSEmin_ind; h2KLmin_ind=fKLmin_ind;
hOffMSEmin_ind=fMSEmin_ind; hOffKLmin_ind=fKLmin_ind;
lMSEhmin_ind=fMSEmin_ind; lKLhmin_ind=fKLmin_ind;lMSEfmin_ind=fMSEmin_ind; lKLfmin_ind=fKLmin_ind;
ldMSEhmin_ind=fMSEmin_ind; ldKLhmin_ind=fKLmin_ind;ldMSEfmin_ind=fMSEmin_ind; ldKLfmin_ind=fKLmin_ind;
fMSEcv=nan(1,subjects); fKLcv= fMSEcv; h1MSEcv=fMSEcv; h1KLcv=fMSEcv; 
h2MSEcv=fMSEcv; h2KLcv=fMSEcv; hOffMSEcv=fMSEcv; hOffKLcv=fMSEcv; 
lMSEfcv=fMSEcv; lMSEhcv=fMSEcv; lKLfcv=fMSEcv; lKLhcv=fMSEcv;
ldMSEfcv=fMSEcv;ldMSEhcv=fMSEcv; ldKLfcv=fMSEcv; ldKLhcv=fMSEcv;
% MSE=cell(1,subjects);KL=cell(1,subjects); 
parpool('local',64)
parfor k=1:subjects

    %preprocess to find relevant subject index (no need for distractor)
    ind=(k-1)*200+1:k*200;
    Tar=wrapToPi(Target(ind) .*(pi/180));
    Dis=wrapToPi(Distractor(ind) .*(pi/180));
    off=wrapToPi(offset(k).*(pi/180));
    Tim=Time(ind);
 
    % mixed bias- uniform for 500 ms trials
    indLong=find(Tim==4000); %short trials
    TarLong=Tar(indLong);
    DisLong=Dis(indLong);
    sResp=wrapToPi(Report(indLong).*(pi/180));
    sDist= (1-cos(sResp-TarLong));
    TrialTrain=round(0.8*length(indLong));
    TrialTest=round(0.2*length(indLong));
    
   
%     for cv=1:5
        disp(['k=' num2str(k) ' cv=' num2str(cv)])
        indTest=cv:5:length(indLong);
        indTrain=1:length(indLong);indTrain(indTest)=[];
        TarTrain=TarLong(indTrain);
        DisTest=DisLong(indTest);
        DisTrain=DisLong(indTrain);
        TarTest=TarLong(indTest);
        sRespTrain=sResp(indTrain);
        sRespTest=sResp(indTest);
        [sHistTrain,~]=histcounts(sRespTrain, edges); sHistTrain=sHistTrain+1; sHistTrain=sHistTrain./(length(indTrain)+length(edges)-1);
        [sHistTest,~]=histcounts(sRespTest, edges); sHistTest=sHistTest+1; sHistTest=sHistTest./(length(indTest)+length(edges)-1);

        %%%% flat potential training
       Resp=nan(TrialTrain,nSim); Dist=nan(TrialTrain,nSim);fMSE=nan(1,length(D_F));fKL=nan(1,length(D_F));
        mses=nan(nSim,1); 
         for dd=1:length(D_F)
             for ns=1:nSim
                 [Resp(:,ns),Dist(:,ns)] = FlatParticle(TrialTrain,D_F(dd),dt,nt,TarTrain);
                 mses(ns)=1/TrialTrain*sum(angle(exp(1i*Resp(:,ns))./exp(1i*sRespTrain)).^2);
             end
             fMSE(dd)=mean(mses);
             [fHist,~]=histcounts(Resp, edges); fHist=fHist+1; fHist=fHist./(TrialTrain*nSim+length(edges)-1);
             fKL(dd)=sum(sHistTrain.*log(sHistTrain./fHist));

         end
         fMSEmin=min(fMSE);
         fMSEmin_ind{k}=find(fMSEmin==fMSE);
         fKLmin=min(fKL);
         fKLmin_ind{k}=find(fKLmin==fKL);

         
         %%%% flat potential testing
       RespM=nan(TrialTest,nSim); DistM=nan(TrialTest,nSim);RespK=nan(TrialTest,nSim); DistK=nan(TrialTest,nSim);
        mses=nan(nSim,1); 
        for ns=1:nSim
                 [RespM(:,ns),DistM(:,ns)] = FlatParticle(TrialTest,D_F(fMSEmin_ind{k}(1)),dt,nt,TarTest);
                 mses(ns)=1/TrialTest*sum(angle(exp(1i*RespM(:,ns))./exp(1i*sRespTest)).^2);
%                  mses(ns)=1/TrialTest*sum((circ_dist(RespM(:,ns),sRespTest)).^2);

                 [RespK(:,ns),DistK(:,ns)] = FlatParticle(TrialTest,D_F(fKLmin_ind{k}),dt,nt,TarTest);
        end
            fMSEcv(k)=mean(mses);
             [fHistK,~]=histcounts(RespK, edges); fHistK=fHistK+1; fHistK=fHistK./(TrialTest*nSim+length(edges)-1);
             fKLcv(k)=sum(sHistTest.*log(sHistTest./fHistK));
      
     
        %%%% het potential 1 train
       Resp=nan(TrialTrain,nSim);Dist=nan(TrialTrain,nSim);
        mses=nan(nSim,1);h1MSE=nan(1,nParam); h1KL=nan(1,nParam);
        het1Param=length(D_H)*length(nTrough1)*length(A1);
        r_ind=randperm(het1Param,nParam);
        for r=1:nParam
            [d,nn,aa]=ind2sub([length(D_H),length(nTrough1),length(A1)],r_ind(r));
            for ns=1:nSim
                [Resp(:,ns),Dist(:,ns)]= HetParticle(TrialTrain,D_H(d),dt,A1(aa),nt,TarTrain,nTrough1(nn), off);
                 mses(ns)=1/TrialTrain*sum(angle(exp(1i*Resp(:,ns))./exp(1i*sRespTrain)).^2);
            end
             h1MSE(r)=mean(mses);
             [h1Hist,~]=histcounts(Resp, edges); h1Hist=h1Hist+1; h1Hist=h1Hist./(TrialTrain*nSim+length(edges)-1);
             h1KL(r)=sum(sHistTrain.*log(sHistTrain./h1Hist));
        end

         h1MSEmin=min(h1MSE);
         h1MSEmin_ind{k}=find(h1MSEmin==h1MSE);
         h1KLmin=min(h1KL);
         h1KLmin_ind{k}=find(h1KLmin==h1KL);

         %%%% het potential 1 test
       RespM=nan(TrialTest,nSim); DistM=nan(TrialTest,nSim);RespK=nan(TrialTest,nSim); DistK=nan(TrialTest,nSim);
       mses=nan(nSim,1); 
        [dM,nnM,aaM]=ind2sub([length(D_H),length(nTrough1),length(A1)],h1MSEmin_ind{k}(1));
        [dK,nnK,aaK]=ind2sub([length(D_H),length(nTrough1),length(A1)],h1KLmin_ind{k}(1));
            for ns=1:nSim
                [RespM(:,ns),DistM(:,ns)]= HetParticle(TrialTest,D_H(dM),dt,A1(aaM),nt,TarTest,nTrough1(nnM), off);
                 mses(ns)=1/TrialTest*sum(angle(exp(1i*RespM(:,ns))./exp(1i*sRespTest)).^2);
                [RespK(:,ns),DistK(:,ns)]= HetParticle(TrialTest,D_H(dK),dt,A1(aaK),nt,TarTest,nTrough1(nnK), off);
            end
             h1MSEcv(k)=mean(mses);
             [h1HistK,~]=histcounts(RespK, edges); h1HistK=h1HistK+1; h1HistK=h1HistK./(TrialTest*nSim+length(edges)-1);
             h1KLcv(k)=sum(sHistTest.*log(sHistTest./h1HistK));
      
       %%%% het potential 2 training
       Resp=nan(TrialTrain,nSim);Dist=nan(TrialTrain,nSim);
        mses=nan(nSim,1);h2MSE=nan(1,nParam); h2KL=nan(1,nParam); 
        het2Param=length(D_H)*length(nTrough1)*length(A1)*length(nTrough2)*length(A2);
        r_ind=randperm(het2Param,nParam);
        for r=1:nParam
            [d,nn1,aa1,nn2,aa2]=ind2sub([length(D_H),length(nTrough1),length(A1),length(nTrough2),length(A2)],r_ind(r));
            for ns=1:nSim
                [Resp(:,ns),Dist(:,ns)]= Het2Particle(TrialTrain,D_H(d),dt,A1(aa1),nt,TarTrain,nTrough1(nn1),A2(aa2),nTrough2(nn2), off);
                 mses(ns)=1/TrialTrain*sum(angle(exp(1i*Resp(:,ns))./exp(1i*sRespTrain)).^2);
            end
             h2MSE(r)=mean(mses);
             [h2Hist,~]=histcounts(Resp, edges); h2Hist=h2Hist+1; h2Hist=h2Hist./(TrialTrain*nSim+length(edges)-1);
             h2KL(r)=sum(sHistTrain.*log(sHistTrain./h2Hist));
        end

         h2MSEmin=min(h2MSE);
         h2MSEmin_ind{k}=find(h2MSEmin==h2MSE);
         h2KLmin=min(h2KL);
         h2KLmin_ind{k}=find(h2KLmin==h2KL);

       %%%% het potential 2 test
       RespM=nan(TrialTest,nSim); DistM=nan(TrialTest,nSim);RespK=nan(TrialTest,nSim); DistK=nan(TrialTest,nSim);
       mses=nan(nSim,1); 
       [dM,nn1M,aa1M,nn2M,aa2M]=ind2sub([length(D_H),length(nTrough1),length(A1),length(nTrough2),length(A2)],h2MSEmin_ind{k}(1));
       [dK,nn1K,aa1K,nn2K,aa2K]=ind2sub([length(D_H),length(nTrough1),length(A1),length(nTrough2),length(A2)],h2KLmin_ind{k}(1));
            for ns=1:nSim
                [RespM(:,ns),DistM(:,ns)]= Het2Particle(TrialTest,D_H(dM),dt,A1(aa1M),nt,TarTest,nTrough1(nn1M),A2(aa2M),nTrough2(nn2M), off);
                 mses(ns)=1/TrialTest*sum(angle(exp(1i*RespM(:,ns))./exp(1i*sRespTest)).^2);
                [RespK(:,ns),DistK(:,ns)]= Het2Particle(TrialTest,D_H(dK),dt,A1(aa1K),nt,TarTest,nTrough1(nn1K),A2(aa2K),nTrough2(nn2K), off);
            end
        h2MSEcv(k)=mean(mses);
        [h2HistK,~]=histcounts(RespK, edges); h2HistK=h2HistK+1; h2HistK=h2HistK./(TrialTest*nSim+length(edges)-1);
        h2KLcv(k)=sum(sHistTest.*log(sHistTest./h2HistK));
            
       %%%% het potential variable offset training
       Resp=nan(TrialTrain,nSim);Dist=nan(TrialTrain,nSim);
       mses=nan(nSim,1);hOffMSE=nan(1,nParam); hOffKL=nan(1,nParam);
       hetOffParam=length(D_H)*length(nTrough1)*length(A1)*length(offs);
       r_ind=randperm(hetOffParam,nParam);
        for r=1:nParam
            [d,nn1,aa1,oo]=ind2sub([length(D_H),length(nTrough1),length(A1),length(offs)],r_ind(r));
            for ns=1:nSim
                [Resp(:,ns),Dist(:,ns)]= HetOffParticle(TrialTrain,D_H(d),dt,A1(aa1),nt,TarTrain,nTrough1(nn1), offs(oo));
                 mses(ns)=1/TrialTrain*sum(angle(exp(1i*Resp(:,ns))./exp(1i*sRespTrain)).^2);
            end
             hOffMSE(r)=mean(mses);
             [hOffHist,~]=histcounts(Resp, edges); hOffHist=hOffHist+1; hOffHist=hOffHist./(TrialTrain*nSim+length(edges)-1);
             hOffKL(r)=sum(sHistTrain.*log(sHistTrain./hOffHist));
        end
         hOffMSEmin=min(hOffMSE);
         hOffMSEmin_ind{k}=find(hOffMSEmin==hOffMSE);
         hOffKLmin=min(hOffKL);
         hOffKLmin_ind{k}=find(hOffKLmin==hOffKL);

        %%%% het potential variable offset testing
         RespM=nan(TrialTest,nSim); DistM=nan(TrialTest,nSim);RespK=nan(TrialTest,nSim); DistK=nan(TrialTest,nSim);
         mses=nan(nSim,1);
         [dM,nn1M,aa1M,ooM]=ind2sub([length(D_H),length(nTrough1),length(A1),length(offs)],hOffMSEmin_ind{k}(1));
         [dK,nn1K,aa1K,ooK]=ind2sub([length(D_H),length(nTrough1),length(A1),length(offs)],hOffKLmin_ind{k}(1));
            for ns=1:nSim
                [RespM(:,ns),DistM(:,ns)]= HetOffParticle(TrialTest,D_H(dM),dt,A1(aa1M),nt,TarTest,nTrough1(nn1M), offs(ooM));
                 mses(ns)=1/TrialTest*sum(angle(exp(1i*RespM(:,ns))./exp(1i*sRespTest)).^2);
                [RespK(:,ns),DistK(:,ns)]= HetOffParticle(TrialTest,D_H(dK),dt,A1(aa1K),nt,TarTest,nTrough1(nn1K), offs(ooK));
            end
        hOffMSEcv(k)=mean(mses);
        [hOffHistK,~]=histcounts(RespK, edges); hOffHistK=hOffHistK+1; hOffHistK=hOffHistK./(TrialTest*nSim+length(edges)-1);
        hOffKLcv(k)=sum(sHistTest.*log(sHistTest./hOffHistK));


        %%%% learning potentials- flat and het training
       l_Resp_H=nan(TrialTrain,nSim);  l_Resp_F=nan(TrialTrain,nSim);
        l_Dist_H=nan(TrialTrain,nSim);l_Dist_F=nan(TrialTrain,nSim);
         msesLH=nan(nSim,1);     msesLF=nan(nSim,1);
         lMSEh=nan(1,nParam); lKLh=nan(1,nParam);lMSEf=nan(1,nParam); lKLf=nan(1,nParam);
        learnParam=length(k_P)*length(s_C)*length(D_L);
        r_ind=randperm(learnParam,nParam);
     for r=1:nParam
 
        [kp,sc,d]=ind2sub([length(k_P),length(s_C),length(D_L)],r_ind(r));
       
        learnPotHet = LearnParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(sc), s_H, k_P(kp));
        
        learnPotFlat = LearnParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(sc), s_H, k_P(kp));

             for ns=1:nSim

               [l_Resp_H(:,ns),l_Dist_H(:,ns)] = LearnParticle(TrialTrain,learnPotHet,D_L(d),dx,dt,xfs,nt,TarTrain, indTrain);
              
               [l_Resp_F(:,ns),l_Dist_F(:,ns)] = LearnParticle(TrialTrain,learnPotFlat,D_L(d),dx,dt,xfs,nt,TarTrain, indTrain);
                msesLH(ns)=1/TrialTrain*sum(angle(exp(1i*l_Resp_H(:,ns))./exp(1i*sRespTrain)).^2);
                msesLF(ns)=1/TrialTrain*sum(angle(exp(1i*l_Resp_F(:,ns))./exp(1i*sRespTrain)).^2);
%                msesLH(ns)=1/TrialTrain*sum((circ_dist(l_Resp_H(:,ns),sRespTrain)).^2);
%                msesLF(ns)=1/TrialTrain*sum((circ_dist(l_Resp_F(:,ns),sRespTrain)).^2);
                
             end         
          lMSEh(r)=mean(msesLH); 
          lMSEf(r)=mean(msesLF); 
          [lHistH,~]=histcounts(l_Resp_H, edges); lHistH=lHistH+1; lHistH=lHistH./(TrialTrain*nSim+length(edges)-1);
          [lHistF,~]=histcounts(l_Resp_F, edges); lHistF=lHistF+1; lHistF=lHistF./(TrialTrain*nSim+length(edges)-1);
          lKLh(r)=sum(sHistTrain.*log(sHistTrain./lHistH));lKLf(r)=sum(sHistTrain.*log(sHistTrain./lHistF));
     end
        lFlatMSEmin=min(lMSEf);
         lMSEfmin_ind{k}=find(lFlatMSEmin==lMSEf);
         lFlatKLmin=min(lKLf);
         lKLfmin_ind{k}=find( lFlatKLmin==lKLf);
          lHetMSEmin=min(lMSEh);
         lMSEhmin_ind{k}=find(lHetMSEmin==lMSEh);
         lHetKLmin=min(lKLh);
         lKLhmin_ind{k}=find( lHetKLmin==lKLh);

         
        %%%% learning flat and het testing
       l_Resp_HM=nan(TrialTest,nSim);  l_Resp_FM=nan(TrialTest,nSim);
        l_Dist_HM=nan(TrialTest,nSim);l_Dist_FM=nan(TrialTest,nSim);
         msesLH=nan(nSim,1);     msesLF=nan(nSim,1);
          l_Resp_HK=nan(TrialTest,nSim);  l_Resp_FK=nan(TrialTest,nSim);
        l_Dist_HK=nan(TrialTest,nSim);l_Dist_FK=nan(TrialTest,nSim);
        [kpFM,scFM,dFM]=ind2sub([length(k_P),length(s_C),length(D_L)],lMSEfmin_ind{k}(1));
        [kpHM,scHM,dHM]=ind2sub([length(k_P),length(s_C),length(D_L)],lMSEhmin_ind{k}(1));
         [kpFK,scFK,dFK]=ind2sub([length(k_P),length(s_C),length(D_L)],lKLfmin_ind{k}(1));
        [kpHK,scHK,dHK]=ind2sub([length(k_P),length(s_C),length(D_L)],lKLhmin_ind{k}(1));
             
        learnPotHetM = LearnParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(scHM), s_H, k_P(kpHM));
        learnPotFlatM = LearnParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(scFM), s_H, k_P(kpFM));
        learnPotHetK = LearnParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(scHK), s_H, k_P(kpHK));
        learnPotFlatK = LearnParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(scFK), s_H, k_P(kpFK));
        tic
            for ns=1:nSim

               [l_Resp_HM(:,ns),l_Dist_HM(:,ns)] = LearnParticle(TrialTest,learnPotHetM,D_L(dHM),dx,dt,xfs,nt,TarTest, indTest);
               [l_Resp_FM(:,ns),l_Dist_FM(:,ns)] = LearnParticle(TrialTest,learnPotFlatM,D_L(dFM),dx,dt,xfs,nt,TarTest, indTest);
               msesLH(ns)=1/TrialTest*sum(angle(exp(1i*l_Resp_HM(:,ns))./exp(1i*sRespTest)).^2);
               msesLF(ns)=1/TrialTest*sum(angle(exp(1i*l_Resp_FM(:,ns))./exp(1i*sRespTest)).^2);
%                msesLH(ns)=1/TrialTest*sum((circ_dist(l_Resp_HM(:,ns),sRespTest)).^2);
%                msesLF(ns)=1/TrialTest*sum((circ_dist(l_Resp_FM(:,ns),sRespTest)).^2);
               [l_Resp_HK(:,ns),l_Dist_HK(:,ns)] = LearnParticle(TrialTest,learnPotHetK,D_L(dHK),dx,dt,xfs,nt,TarTest, indTest);
               [l_Resp_FK(:,ns),l_Dist_FK(:,ns)] = LearnParticle(TrialTest,learnPotFlatK,D_L(dFK),dx,dt,xfs,nt,TarTest, indTest);

            end
            toc
            
        lMSEfcv(k)=mean(msesLF);
        [lHistKf,~]=histcounts(l_Resp_FK, edges); lHistKf=lHistKf+1; lHistKf=lHistKf./(TrialTest*nSim+length(edges)-1);
        lKLfcv(k)=sum(sHistTest.*log(sHistTest./lHistKf));
        lMSEhcv(k)=mean(msesLH);
        [lHistKh,~]=histcounts(l_Resp_HK, edges); lHistKh=lHistKh+1; lHistKh=lHistKh./(TrialTest*nSim+length(edges)-1);
        lKLhcv(k)=sum(sHistTest.*log(sHistTest./lHistKh));
            
        
          %%%% learning distractor potentials- flat and het training
       l_Resp_H=nan(TrialTrain,nSim);  l_Resp_F=nan(TrialTrain,nSim);
        l_Dist_H=nan(TrialTrain,nSim);l_Dist_F=nan(TrialTrain,nSim);
         msesLH=nan(nSim,1);     msesLF=nan(nSim,1);
         lMSEh=nan(1,nParam); lKLh=nan(1,nParam);lMSEf=nan(1,nParam); lKLf=nan(1,nParam);
        learnParam=length(k_P)*length(s_C)*length(D_L);
        r_ind=randperm(learnParam,nParam);
     for r=1:nParam
         tic
        [kp,sc,d]=ind2sub([length(k_P),length(s_C),length(D_L)],r_ind(r));
        learnDisPotHet = LearnDistParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(sc), s_H, k_P(kp), DisLong);
        learnDisPotFlat = LearnDistParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(sc), s_H, k_P(kp), DisLong);
        toc
        tic
             for ns=1:nSim
               [l_Resp_H(:,ns),l_Dist_H(:,ns)] = LearnDistParticle(TrialTrain,D_L(d),dx,dt,xfs,nt,TarTrain,learnDisPotHet,DisTrain, indTrain);
               [l_Resp_F(:,ns),l_Dist_F(:,ns)] = LearnDistParticle(TrialTrain,D_L(d),dx,dt,xfs,nt,TarTrain,learnDisPotFlat,DisTrain, indTrain);
                msesLH(ns)=1/TrialTrain*sum(angle(exp(1i*l_Resp_H(:,ns))./exp(1i*sRespTrain)).^2);
                msesLF(ns)=1/TrialTrain*sum(angle(exp(1i*l_Resp_F(:,ns))./exp(1i*sRespTrain)).^2);
             end
             toc
          lMSEh(r)=mean(msesLH); 
          lMSEf(r)=mean(msesLF); 
          [lHistH,~]=histcounts(l_Resp_H, edges); lHistH=lHistH+1; lHistH=lHistH./(TrialTrain*nSim+length(edges)-1);
          [lHistF,~]=histcounts(l_Resp_F, edges); lHistF=lHistF+1; lHistF=lHistF./(TrialTrain*nSim+length(edges)-1);
          lKLh(r)=sum(sHistTrain.*log(sHistTrain./lHistH));lKLf(r)=sum(sHistTrain.*log(sHistTrain./lHistF));
     end
        ldFlatMSEmin=min(lMSEf);
         ldMSEfmin_ind{k}=find(ldFlatMSEmin==lMSEf);
         ldFlatKLmin=min(lKLf);
         ldKLfmin_ind{k}=find( ldFlatKLmin==lKLf);
          ldHetMSEmin=min(lMSEh);
         ldMSEhmin_ind{k}=find(ldHetMSEmin==lMSEh);
         ldHetKLmin=min(lKLh);
         ldKLhmin_ind{k}=find( ldHetKLmin==lKLh);
         
            
        %%%% learning distractor flat and het testing
       l_Resp_HM=nan(TrialTest,nSim);  l_Resp_FM=nan(TrialTest,nSim);
        l_Dist_HM=nan(TrialTest,nSim);l_Dist_FM=nan(TrialTest,nSim);
         msesLH=nan(nSim,1);     msesLF=nan(nSim,1);
          l_Resp_HK=nan(TrialTest,nSim);  l_Resp_FK=nan(TrialTest,nSim);
        l_Dist_HK=nan(TrialTest,nSim);l_Dist_FK=nan(TrialTest,nSim);
        [kpFM,scFM,dFM]=ind2sub([length(k_P),length(s_C),length(D_L)],ldMSEfmin_ind{k}(1));
        [kpHM,scHM,dHM]=ind2sub([length(k_P),length(s_C),length(D_L)],ldMSEhmin_ind{k}(1));
         [kpFK,scFK,dFK]=ind2sub([length(k_P),length(s_C),length(D_L)],ldKLfmin_ind{k}(1));
        [kpHK,scHK,dHK]=ind2sub([length(k_P),length(s_C),length(D_L)],ldKLhmin_ind{k}(1));
             
        learnPotHetM = LearnDistParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(scHM), s_H, k_P(kpHM), DisLong);
        learnPotFlatM = LearnDistParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(scFM), s_H, k_P(kpFM), DisLong);
        learnPotHetK = LearnDistParticleLandscape(Trials,baselinePot_het,dx,xfs,TarLong, s_C(scHK), s_H, k_P(kpHK), DisLong);
        learnPotFlatK = LearnDistParticleLandscape(Trials,baselinePot_flat,dx,xfs,TarLong, s_C(scFK), s_H, k_P(kpFK), DisLong);
            for ns=1:nSim
               [l_Resp_HM(:,ns),l_Dist_HM(:,ns)] = LearnDistParticle(TrialTest,D_L(dHM),dx,dt,xfs,nt,TarTest,learnPotHetM, DisTest, indTest);
               [l_Resp_FM(:,ns),l_Dist_FM(:,ns)] = LearnDistParticle(TrialTest,D_L(dFM),dx,dt,xfs,nt,TarTest, learnPotFlatM, DisTest,indTest);
               msesLH(ns)=1/TrialTest*sum(angle(exp(1i*l_Resp_HM(:,ns))./exp(1i*sRespTest)).^2);
               msesLF(ns)=1/TrialTest*sum(angle(exp(1i*l_Resp_FM(:,ns))./exp(1i*sRespTest)).^2);
               [l_Resp_HK(:,ns),l_Dist_HK(:,ns)] = LearnDistParticle(TrialTest,D_L(dHK),dx,dt,xfs,nt,TarTest, learnPotHetK, DisTest, indTest);
               [l_Resp_FK(:,ns),l_Dist_FK(:,ns)] = LearnDistParticle(TrialTest,D_L(dFK),dx,dt,xfs,nt,TarTest,learnPotFlatK, DisTest, indTest);
            end
            
        ldMSEfcv(k)=mean(msesLF);
        [lHistKf,~]=histcounts(l_Resp_FK, edges); lHistKf=lHistKf+1; lHistKf=lHistKf./(TrialTest*nSim+length(edges)-1);
        ldKLfcv(k)=sum(sHistTest.*log(sHistTest./lHistKf));
        ldMSEhcv(k)=mean(msesLH);
        [lHistKh,~]=histcounts(l_Resp_HK, edges); lHistKh=lHistKh+1; lHistKh=lHistKh./(TrialTest*nSim+length(edges)-1);
        ldKLhcv(k)=sum(sHistTest.*log(sHistTest./lHistKh));

    
    
%     MSE{k}(1)=mean(fMSEcv);
%     KL{k}(1)=mean(fKLcv);
%     MSE{k}(2)=mean(h1MSEcv);
%     KL{k}(2)=mean(h1KLcv);
%     MSE{k}(3)=mean(h2MSEcv);
%     KL{k}(3)=mean(h2KLcv);
%     MSE{k}(4)=mean(hOffMSEcv);
%     KL{k}(4)=mean(hOffKLcv);
%     MSE{k}(5)=mean(lMSEfcv);
%     KL{k}(5)=mean(lKLfcv);
%     MSE{k}(6)=mean(lMSEhcv);
%     KL{k}(6)=mean(lKLhcv);
%     MSE{k}(7)=mean(ldMSEfcv);
%     KL{k}(7)=mean(ldKLfcv);
%     MSE{k}(8)=mean(ldMSEhcv);
%     KL{k}(8)=mean(ldKLhcv);


end

delete(gcp('nocreate'))


save('LongTrialModelFits_CV1','fKLmin_ind','fMSEmin_ind','h1KLmin_ind','h1MSEmin_ind','h2KLmin_ind','h2MSEmin_ind',...
    'hOffKLmin_ind','hOffMSEmin_ind','lKLfmin_ind','lMSEfmin_ind','lKLhmin_ind','lMSEhmin_ind',...
   'ldKLfmin_ind','ldMSEfmin_ind','ldKLhmin_ind','ldMSEhmin_ind','fMSEcv','fKLcv', 'h1MSEcv', 'h1KLcv',...
'h2MSEcv', 'h2KLcv', 'hOffMSEcv', 'hOffKLcv', 'lMSEfcv', 'lMSEhcv', 'lKLfcv', 'lKLhcv',...
'ldMSEfcv', 'ldMSEhcv', 'ldKLfcv', 'ldKLhcv','-v7.3')

