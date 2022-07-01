function [fResp,fDist] = FlatParticle(Trials,D_F,dt,nt,Tar )
%computes the diffusion of a flat potential landscape for given trials
%  %Trials= number of trials,  Df=
%diffusion, dx=response stepsize, dt= timestep,  xfs= values of response, nt= number
%of timesteps, Tar= targets, off= het offset,

%initialize responses and distortion
    fResp=zeros(Trials,1);
    fDist=fResp;
   
   for tr=1:Trials
       inp=Tar(tr); %target
       
       %initialize particle
       flatPart=inp;

       %particle trajectory
       for j=1:nt-1
           flatPart=mod(flatPart+sqrt(dt*2*D_F)*randn+pi,2*pi)-pi;

       end
         
     %model responses and distortion for all
     fResp(tr)=flatPart;
     fDist(tr)=(1-cos(fResp(tr)-inp));
   end
end

