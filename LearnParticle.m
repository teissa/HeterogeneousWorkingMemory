function [lResp,lDist] = LearnParticle(Trials,learnPot,D_L,dx,dt,xfs,nt,Tar,ind)
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
    lResp=zeros(Trials,1);
    lDist=lResp;

   for tr=1:Trials
       inp=Tar(tr); %target

       
       %initialize particles
       learnPart=inp;

       
       %gradient of learning model
%        learnGrad=gradient(learnPot(ind(tr),:),dx);


learnGrad=diff(learnPot(ind(tr),:))./dx;
learnGrad=[learnGrad,learnGrad(1)];


       %particle trajectory

       for j=1:nt-1
           learnPart=mod(learnPart-dt*interp1(xfs,learnGrad,learnPart)+sqrt(dt*2*D_L)*randn+pi,2*pi)-pi;

       end

         
         %model responses and distortion for all

         lResp(tr)=learnPart;
         lDist(tr)=(1-cos(lResp(tr)-inp));
         clear learnGrad
   end
   
end

