clear
close all

D = 0;   % noise corresponding to corrosion of parametric WM rep, diffusion
TT = [0.1,0.5,1,5,10];      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
pr = @(x) exp(Am*cos(m*x)); %prior function
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

xfs = linspace(-pi,pi,256); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size


% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
ns=1:12; %number of troughs in the potential 
As=1; %As = linspace(0,2,21);   % varying amplitudes of the heterogeneous potential 

close all
for t=1:length(TT)
    T=TT(t);
    nt = round(T/dt)+1; % number of timesteps per sim
  rdist{t} = zeros(length(ns),length(As)); %initialize distortion
for r=1:length(ns)
    n=ns(r)
for k=1:length(As)
    A = As(k); %amplitude of potential for this iteration
    disp(['Working on A=' num2str(A) ' n=' num2str(n)])
    for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=inp; %starting point of particle
        for j=1:nt-1
            x=mod(x-dt*A*sin(n*x)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        rdist{t}(r,k) = rdist{t}(r,k)+(1-cos(x-inp))/Nsim;  % running average of distortion
    end
end
end
end

%load diffDist_varDT

figure; hold on;
clrs=[linspace(0,0.8,5);linspace(0,0.8,5);linspace(0,0.8,5)]';
for t=1:5
    plot(rdist{t}(1:12),'linewidth',5,'color',clrs(t,:))
end
set(gca,'fontsize',30);
set(gca,'ticklabelinterpreter','latex')
xlabel('Wells','interpreter','latex');
ylabel('Mean Distortion','interpreter','latex')
legend('0.1s','0.5s','1s','5s','10s','interpreter','latex')

%%

clear all
close all

DD = [0.01,0.05,0.1];   % noise corresponding to corrosion of parametric WM rep, diffusion
TT = [0.1,0.5,1,5,10];      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
pr = @(x) exp(Am*cos(m*x)); %prior function
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function
Nsim=1e5;
randinp = zeros(1,Nsim); % this identifies the angle theta for each simulation


for d=1:length(DD)
    D=DD(d);
    for t=1:length(TT)
        T=TT(t);
        % average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
ns=1:12; %number of troughs in the potential 
As = 1;   % varying amplitudes of the heterogeneous potential 
rdist{d,t} = zeros(length(ns),length(As)); %initialize distortion
clrs=parula(length(ns));

        xfs = linspace(-pi,pi,201); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size

for r=1:length(ns)
    n=ns(r);
for k=1:length(As)
    A = As(k); %amplitude of potential for this iteration
    disp(['Working on A=' num2str(A) ' n=' num2str(n)])
    for l=1:Nsim
        inp = randinp(l); %angle theta from distribution creted above
        x=inp; %starting point of particle
        for j=1:nt-1
            x=mod(x-dt*A*sin(n*x)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
        rdist{d,t}(r,k) = rdist{d,t}(r,k)+(1-cos(x-inp))/Nsim;  % running average of distortion
    end
end
end
    end
end

%load diffDist_varDT

figure; hold on;
clrs=[linspace(0,0.8,5);linspace(0,0.8,5);linspace(0,0.8,5)]';
for t=1:5
    plot(rdist{2,t},'linewidth',5,'color',clrs(t,:))
end
set(gca,'fontsize',30);
set(gca,'ticklabelinterpreter','latex')
xlabel('Wells','interpreter','latex');
ylabel('Mean Distortion','interpreter','latex')
legend('0.1s','0.5s','1s','5s','10s','interpreter','latex')

%%%%% SFig


figure; hold on;
clrs=[linspace(0,0.8,5);linspace(0,0.8,5);linspace(0,0.8,5)]';
for t=1:5
    plot(rdist{1,t},'linewidth',5,'color',clrs(t,:))
end
set(gca,'fontsize',30);
set(gca,'ticklabelinterpreter','latex')
xlabel('Wells','interpreter','latex');
ylabel('Mean Distortion','interpreter','latex')
legend('0.1s','0.5s','1s','5s','10s','interpreter','latex')

figure; hold on;
clrs=[linspace(0,0.8,5);linspace(0,0.8,5);linspace(0,0.8,5)]';
for t=1:5
    plot(rdist{3,t},'linewidth',5,'color',clrs(t,:))
end
set(gca,'fontsize',30);
set(gca,'ticklabelinterpreter','latex')
xlabel('Wells','interpreter','latex');
ylabel('Mean Distortion','interpreter','latex')
legend('0.1s','0.5s','1s','5s','10s','interpreter','latex')
