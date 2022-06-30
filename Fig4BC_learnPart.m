%use trials to learn the heterogeneous potential landscape for a particle
%model. This is a simplified version of what we will use in future
%iterations
clear
close all

%create the prior and parameterize it to build the sampling 

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
pr = @(x) exp(Am*cos(m*x)); %prior function
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

% build sampling
Nsim = 500;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 5e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from 
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation

kap =8;  % amplitude of von mises to be used
sh = 1/4;   % shift
sc = 5;     % scaling
xs = linspace(-pi,pi,1000); dx = xs(2)-xs(1);

% flat potential
fsol = 1/2/pi;
for n=1:10, fsol = fsol+cos(n*xs).*exp(-D*n^2*T)/pi; end
fdist = dx*sum((1-cos(xs)).*fsol);


%initialize parameters
pot = 1/2/pi+0*xs; %initial potential landscape (flat)
 figure, hold on,  xlim([-180 180]); yticks([])
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');
plot(xs.*(180/pi),pot,'linewidth',4,'Color','k');
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')

for k=1:10
pot(k,:) = 1/2/pi+0*xs; %initial potential landscape (flat)
% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation

rout = zeros(Nsim,1); %responses at end of trial
rdist = 0; %distortion
rdistHet=0;
n=m; A=1; %het model
for l=1:Nsim
    inp = randinp(l);
    x=inp; xhet=inp;
    pvar(l,:)=pot(k,:);

    pgrad = gradient(pot(k,:),dx); %find the derivative of the current potential landscape
    
    % simulate the delayed response trial
    for j=1:nt-1 
        x=mod(x-dt*interp1(xs,pgrad,x)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; 
        xhet=mod(xhet-dt*A*sin(n*xhet)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; 
    end
    rout(l) = x;%response at end with drift and diffusion
    
    % add the distortion value for this trial
    rdist = (rdist+(1-cos(x-inp)));  % running average of distortion
    rdistHet = (rdistHet+(1-cos(xhet-inp)));  % running average of distortion
    trial_dist(l)=(1-cos(x-inp));
    mdist(k,l)=rdist/l;
    mdistHet(k,l)=rdistHet/l;

    %update the potential landscape using the previous trial input
    pot(k,:) = pot(k,:)+sc*(sh-exp(kap*(cos(xs-inp)-1)))/l;
    pot(k,:) = pot(k,:)/dx/sum(pot(k,:));
    
end
end

%%

figure; hold on,  plot([1:Nsim],0*[1:Nsim]+fdist,'k--','linewidth',2);
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex')
 xlim([1 Nsim]); ylim([0 0.5])
 for k=1:10
plot(mdist(k,:),'Color',[0.5,0.5,0.5],'LineWidth',2); 
 end
 plot(mean(mdist),'k','LineWidth',5); 
plot([1 Nsim],[mean(mean(mdistHet)) mean(mean(mdistHet))],'k-.','LineWidth',2);
xlabel('Trial','interpreter','latex','fontsize',30);
ylabel('Distortion','interpreter','latex','fontsize',30);

%%

figure;  pcolor(xs,[1:Nsim],pvar);  shading flat; colorbar
 xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('Trial','fontsize',30,'interpreter','latex'); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex'); set(gca,'Ydir','reverse')
%% 

 figure, hold on,  xlim([-180 180]); yticks([])
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$U(\theta|\theta_{1:500})$','fontsize',30,'interpreter','latex');
for k=1:10
plot(xs.*(180/pi),pot(k,:),'linewidth',2,'Color',[0.5,0.5,0.5]);
end
plot(xs.*(180/pi),mean(pot),'linewidth',4,'Color','k');
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')

