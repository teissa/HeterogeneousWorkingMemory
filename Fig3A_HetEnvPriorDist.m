%% Flat Landscape Distortion

clear
close all


D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior and landscape
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
n=0;
A=0;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A*cos(n*x);
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
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance

        end
                dist(l)=1-cos(x(end)-inp);

end

figure; plot(xfs, pr(xfs)/sum(pr(xfs)),'r','linewidth',4)
xlabel('$\theta$','Interpreter','Latex')
ylabel('$P_{env}(\theta)$','Interpreter','Latex'); %yticks([]);
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})

figure; plot(xfs,ut(xfs),'b','linewidth',4)
xlabel('$\theta$','Interpreter','Latex')
ylabel('$U(\theta)$','Interpreter','Latex'); yticks([]);
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})

edges=linspace(-pi,pi,30);
xboot=1000;
figure; hold on;
for k=1:29
    ind=find(randinp>=edges(k) & randinp<edges(k+1));
    for j=1:xboot
        r=randi(length(ind),[1,length(ind)]);
        mdistboot(j)=mean(dist(ind(r)));
    end

     mdist(k)=mean(dist(ind));
    sdist(k)=std(mdistboot);
end
bar(edges(1:end-1),mdist,'k')
errorbar(edges(1:end-1),mdist,sdist,'k.')
ylabel('$\bar{d}(\theta)$','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})
ylim([0 0.5])



%% Environment-Matched 4-well landscpae
close all

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
n=4;
A=1;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x);
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
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1))+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
                dist(l)=1-cos(x(end)-inp);

end

figure; hold on; plot(xfs,ut(xfs),'b','linewidth',4)
xlabel('$\theta$','Interpreter','Latex')
ylabel('$U(\theta)$','Interpreter','Latex'); yticks([]); plot([-pi/2 -pi/2],[-1 1],'k--');plot([0 0],[-1 1],'k--');plot([pi/2 pi/2],[-1 1],'k--');
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})



edges=linspace(-pi,pi,30);
xboot=1000;
figure; hold on;
for k=1:29
    ind=find(randinp>=edges(k) & randinp<edges(k+1));
    for j=1:xboot
        r=randi(length(ind),[1,length(ind)]);
        mdistboot(j)=mean(dist(ind(r)));
    end
    mdist(k)=mean(mdistboot);
    sdist(k)=std(mdistboot);
end
bar(edges(2:end),mdist,'k')
errorbar(edges(2:end),mdist,sdist,'k.')
ylabel('$\bar{d}(\theta)$','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})
ylim([0 0.5])

