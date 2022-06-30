% note, spatial units are in radians and time units are in seconds. noise
% amplitude is set to roughly match psychophysics results from various
% papers (Bayes, Wimmer, etc)
clear all
close all

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
n=0;
A=1;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x);
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

Nsim =1e5;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randind=randi(256,1,Nsim);
randinpOn = zeros(1,Nsim); %xsamp(randind); % this identifies the angle theta for each simulation
randinpOff = repmat(pi/4,[1,Nsim]); %xsamp(randind); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,256); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size
% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;
cmap=colormap('hsv');
cmap=cmap([226:256,1:225],:);
colormap(cmap)
for l=1:Nsim
        inp = randinpOn(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        respOn(l)=x(end);
        distOn(l)=1-cos(x(end)-inp);
        
        inp = randinpOff(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance

        end
        respOff(l)=x(end);
        distOff(l)=1-cos(x(end)-inp);
       
end
mdistOff=mean(distOff)
mdistOn=mean(distOn)
edges=linspace(-pi,pi,30);
figure; hold on;
xlabel('Target','Interpreter','Latex')
ylabel('Response Frequency','Interpreter','Latex')
title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); ylim([0 0.5]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
histogram(respOn,edges,'Normalization','probability','FaceColor',[0.75, 0.75,0.75]);%pdOn=fitdist(respOn','Normal'); gaussOn=normpdf(xfs,pdOn.mu,pdOn.sigma); plot(xfs,gaussOn);
histogram(respOff,edges,'Normalization','probability','FaceColor',[0.5, 0.5,0.5]);
plot([0,0],[0,0.5],'LineWidth',5, 'Color',[0.75, 0.75,0.75]);
plot([pi/4,pi/4],[0,0.5],'LineWidth',5, 'Color',[0.5, 0.5,0.5]);

 %plot(randinp,dist,'k.','MarkerSize',20)

%%

n=4;
A=0.4;
Nsim =1e5;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randind=randi(256,1,Nsim);
randinpOn = zeros(1,Nsim); %xsamp(randind); % this identifies the angle theta for each simulation
randinpOff = repmat(pi/4,[1,Nsim]); %xsamp(randind); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,256); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size
% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;
cmap=colormap('hsv');
cmap=cmap([226:256,1:225],:);
colormap(cmap)
for l=1:Nsim
        inp = randinpOn(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
           
        end
        respOn(l)=x(end);
        distOn(l)=1-cos(x(end)-inp);
        
        inp = randinpOff(l); %angle theta from distribution creted above
        x=zeros(1,nt); %starting point of particle
        x(1)=inp;
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        respOff(l)=x(end);
        distOff(l)=1-cos(x(end)-inp);
       
end
mdistOff=mean(distOff)
mdistOn=mean(distOn)
edges=linspace(-pi,pi,30);
figure; hold on;
xlabel('Target','Interpreter','Latex')
ylabel('Response Frequency','Interpreter','Latex')
title('Heterogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); ylim([0 0.5]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
histogram(respOn,edges,'Normalization','probability','FaceColor',[0.75, 0.75,0.75]);%pdOn=fitdist(respOn','Normal'); gaussOn=normpdf(xfs,pdOn.mu,pdOn.sigma); plot(xfs,gaussOn);
histogram(respOff,edges,'Normalization','probability','FaceColor',[0.5, 0.5,0.5]);
plot([0,0],[0,0.5],'LineWidth',5, 'Color',[0.75, 0.75,0.75]);
plot([pi/4,pi/4],[0,0.5],'LineWidth',5, 'Color',[0.5, 0.5,0.5]);
plot([-pi -pi],[0 0.5],'k--');plot([pi pi],[0 0.5],'k--');
plot([-pi/2 -pi/2],[0 0.5],'k--');plot([pi/2 pi/2],[0 0.5],'k--');
plot([0 0],[0 0.5],'k--');