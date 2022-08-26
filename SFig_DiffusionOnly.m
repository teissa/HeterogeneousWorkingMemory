%% Code to analyze just the diffusion effects

clear all
close all

DD = [0.01,0.05,0.1];   % noise corresponding to corrosion of parametric WM rep, diffusion
TT = [0.1,0.5,1,5,10];      % delay time (seconds)

Nsim=1e5;


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

        xfs = linspace(-pi,pi,360); %angles discretized to find distortion
dx = xfs(2)-xfs(1); %angle step size

for r=1:length(ns)
    n=ns(r);

for k=1:length(As)
    A = As(k); %amplitude of potential for this iteration
        ut=@(x)-A*cos(n*x);
        %[mx, mind]=min(ut(xfs));
        mind=2*pi*(round(0.5*n)/n)-pi;
    randinp = repmat(mind,[1,Nsim]); % this identifies the angle theta for each simulation

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