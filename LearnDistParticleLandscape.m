function [learnPot] = LearnDistParticleLandscape(Trials,learnPotBase,dx,xfs,Tar, sc, sh, kap, Dis)
%this function runs the larticle models for a het, flat, and learning
%particle for just the target points of Panichello 19 exp 2. It reports
%back the responses and ditortion for the subject, and all 3 models, as
%well as the learned potential landscape.

%Trials= number of trials, learnPot= initialized learning potential, D*=
%diffusion, dx=response stepsize, dt= timestep,A=het amplitude,  xfs= values of response, nt= number
%of timesteps, Tar= targets, nTrough=troughs for het model, sc/sh/kap=
%parameters for learning, Resp= subject responses, off= het offset,
%subPot=initialized subject potential, learnDistPot= distractor included,
%distractor= distractor targets, fg= figure toggle
  



%initialize responses and distortion
learnPot=nan(Trials,length(learnPotBase));
learnPot(1,:)=learnPotBase;
   for tr=2:Trials
       inp=Tar(tr-1); %target
       dis=Dis(tr-1);
        %update learning model
         learnPot(tr,:) = learnPot(tr-1,:)+(sc*(sh-exp(kap*(cos(xfs-inp)-1)))/tr);
         learnPot(tr,:) = learnPot(tr,:)/dx/sum(learnPot(tr,:));
         learnPot(tr,:) = learnPot(tr-1,:)+(sc*(sh-exp(kap*(cos(xfs-dis)-1)))/tr);
         learnPot(tr,:) = learnPot(tr,:)/dx/sum(learnPot(tr,:));
         

   end
   
end
