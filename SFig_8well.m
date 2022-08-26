%% 8-well offset landscape and particle trajectories


D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
n=8;
A=1;
offset=2.5;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x-offset);
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

Nsim = 100;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,256); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size
% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;
figure; hold on;
xlabel('Time(s)','Interpreter','Latex')
ylabel('Target($\theta$)','Interpreter','Latex')
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
axis([ 0 T -pi pi])
for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1)+offset)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
            while abs(x(j)-x(j-1)) >3 
                 x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1)+offset)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
            end
        end
                dist(l)=1-cos(x(end)-inp);

        plot(tim,x,'Color','k')
end
prb=pr(xfs); prb=(prb-min(prb))./(max(prb)-min(prb));
plot(prb,xfs,'r','linewidth',2)
yticks([-2,0,2]); yticklabels({'-100','0','100'})

U=-ut(xfs);
figure; hold on; plot(xfs,ut(xfs),'k','linewidth',4);
plot([-pi/2 -pi/2],[-1 1],'k--');plot([0 0],[-1 1],'k--');plot([pi/2 pi/2],[-1 1],'k--');
xlabel('$\theta$','Interpreter','Latex')
ylabel('$U(\theta)$','Interpreter','Latex'); yticks([]);
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})

%% distortion


D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
n=8;
A=1;
offset=pi/4;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x-offset);
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

Nsim = 1e5;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,sampres); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size
% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;
for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1)+offset)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
                dist(l)=1-cos(x(end)-inp);

end

edges=linspace(-pi,pi,20);
xboot=1000;
figure; hold on;
for k=1:19
    ind=find(randinp>=edges(k) & randinp<edges(k+1));
    for j=1:xboot
        r=randi(length(ind),[1,length(ind)]);
        mdistboot(j)=mean(dist(ind(r)));
    end
    mdist(k)=mean(mdistboot);
    sdist(k)=std(mdistboot);
end
bar(mdist,'k')
errorbar(mdist,sdist,'k.')
ylabel('Distortion','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
 xlim([1,19]);xticks([])
ylim([0 0.5])

%% conditional and marginal distributions

load P_resp_theta_HetEightOff
figure; hold on;
pcolor(xfs(1:end-1),xfs,condist); shading flat
axis('tight'); colormap('hot')
title('$P(\theta_{resp}|\theta_{env})$','interpreter','latex')
xlabel('$\theta_{resp}$','Interpreter','Latex')
ylabel('$\theta_{env}$','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
set(gca,'box','off'); caxis([0 0.015]); colorbar
yticks([-2,0,2]); yticklabels({'-100','0','100'})
xticks([-2,0,2]); xticklabels({'-100','0','100'})


margdist=sum(condist)/(length(xfs)-1).*pr(xfs(1:end-1));
mardist=margdist/sum(margdist);

figure; hold on; plot(xfs(1:end-1),mardist,'k','linewidth',4)
xlabel('$\theta$','Interpreter','Latex')
ylabel('$P(\theta_{resp})$','Interpreter','Latex'); 
ylim([0 0.025])
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})

%% Builds distortion heatmaps 
close all

close all
count=1;

for tt=1:length(Ts)
        T=Ts(tt);
   for dd=1:length(Ds)
    D=Ds(dd);

% find distortion with flat poential (no troughs)

xfs = linspace(-pi,pi,201); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size

% test a heterogenous potential landscape

% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
ns=1:16; %number of troughs in the potential 
As = linspace(0,2,21);   % varying amplitudes of the heterogeneous potential 
rdist{dd,tt} = zeros(length(ns),length(As)); %initialize distortion

for r=1:length(ns)
    n=ns(r);
for k=1:length(As)
    A = As(k); %amplitude of potential for this iteration
    disp(['Working on A=' num2str(A) ' n=' num2str(n) ' d=' num2str(D) ' T=' num2str(T)])
    for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=inp; %starting point of particle
        for j=1:nt-1
            x=mod(x-dt*A*sin(n*x-2.5)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        rdist{dd,tt}(r,k) = rdist{dd,tt}(r,k)+(1-cos(x-inp))/Nsim;  % running average of distortion
    end
end
end
 count=count+1;
    end
end





%%  Just makes the figures
close all

load SFigOffset_varDT
rOffset=rdist;
load hetP_varDT
d=2;
for k=1:3
figure; hold on; 
pcolor(As,ns,rOffset{d,k}-rdist{d,k}); shading flat; colormap(flipud(hot)); caxis([-0.05 0.05]); colorbar
 [~,mind]=min(rOffset{d,k}-rdist{d,k});
% plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Wells','fontsize',30,'interpreter','latex');
end