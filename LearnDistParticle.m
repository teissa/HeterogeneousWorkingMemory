function [ldResp,ldDist] = LearnDistParticle(Trials,D_L,dx,dt,xfs,nt,Tar,learnDistPot, Distractor,ind)
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
    ldResp=zeros(Trials,1);
    ldDist=ldResp;

   for tr=1:Trials
       inp=Tar(tr); %target
       dis=Distractor(tr);
       
       %initialize particles
       learnDistPart=inp;
       
       %gradient of learning model
      % learnDistGrad=gradient(learnDistPot(ind(tr),:),dx);
      
learnDistGrad=diff(learnDistPot(ind(tr),:))./dx;
learnDistGrad=[learnDistGrad,learnDistGrad(1)];
       %particle trajectory
       
       for j=1:nt-1
           learnDistPart=mod(learnDistPart-dt*interp1(xfs,learnDistGrad,learnDistPart)+sqrt(dt*2*D_L)*randn+pi,2*pi)-pi;

       end

         
         %model responses and distortion for all
         ldResp(tr)=learnDistPart;
         ldDist(tr)=(1-cos(ldResp(tr)-inp));
         clear learnDistGrad
   end
   
end

