function [hResp,hDist] = HetOffParticle(Trials,D_H,dt,A,nt,Tar,nTrough, off)
%computes the diffusion of a heterogeneous potential landscape for given trials
%Trials= number of trials, learnPot= initialized learning potential, D*=
%diffusion, dx=response stepsize, dt= timestep,A=het amplitude,  xfs= values of response, nt= number
%of timesteps, Tar= targets, nTrough=troughs for het model, sc/sh/kap=
%parameters for learning, Resp= subject responses, off= het offset,
%subPot=initialized subject potential, learnDistPot= distractor included,
%distractor= distractor targets, fg= figure toggle
  



%initialize responses and distortion
    hResp=zeros(Trials,1); 
    hDist=hResp;
   
   for tr=1:Trials
       inp=Tar(tr); %target
       
       %initialize particles
       hetPart=inp; 
       
       %particle trajectory
       for j=1:nt-1
           hetPart=mod(hetPart-dt*A*sin(nTrough*hetPart-off)+sqrt(dt*2*D_H)*randn+pi,2*pi)-pi;     
       end
         
         %model responses and distortion for all
         hResp(tr)=hetPart;
         hDist(tr)=(1-cos(hResp(tr)-inp));
   end
   
end

