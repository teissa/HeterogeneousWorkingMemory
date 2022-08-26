% runs sims for conditional prob

clear all
close all

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 1;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
n=4;
A=1;
offset=0;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x-offset);
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)

Nsim = 1e5;   % number of sims per amplitude
% average distortion in each case
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim
tim=0:dt:T;


xfs = linspace(-pi,pi,256); %angles discretized to find distortion
condist=zeros(length(xfs),length(xfs)-1);
for k=1:length(xfs)
    disp(k)
    for s=1:Nsim
        x=zeros(1,nt); %starting point of particle
        x(1)=xfs(k);
        for j=2:nt
            x(j)=mod(x(j-1)-dt*A*sin(n*x(j-1)+offset)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; %movement of the particle, normalized by circular distance
        end
          resps(s)=x(end);
    end
condist(k,:)=histcounts(resps,xfs)./Nsim;
end

%% Figures for Flat Landscape
close all

load P_resp_theta_FLAT
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
xlabel(' Target ($\theta$)','Interpreter','Latex')
ylabel('$P(\theta_{resp})$','Interpreter','Latex'); 
ylim([0 0.025])
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})

%% Figure for Environment-Matched Landscape


load P_resp_theta_HetFour
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
xlabel('Target ($\theta$)','Interpreter','Latex')
ylabel('$P(\theta_{resp})$','Interpreter','Latex'); 
ylim([0 0.025])
%title('Homogeneous','Interpreter','Latex')
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
xlim([-pi pi]); set(gca,'box','off'); xticks([-2,0,2]); xticklabels({'-100','0','100'})




